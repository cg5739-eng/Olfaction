"""2-D passive scalar transport and four-inlet robot. Run: python odor_sim.py"""
import argparse
import numpy as np
from scipy.fft import dctn, idctn
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, RadioButtons, Slider
from matplotlib.patches import Circle, FancyArrowPatch

L=5.; N=100; DX=L/N; D=0.005; EMISSION=1.; DECAY=0.002
DT=0.025; MAX_SPEED=0.2; INLET=0.15; REACH=0.25
BODY_RADIUS=0.15; CAPTURE_RADIUS=0.20; SUCTION=0.012; PORT_SIGMA=0.055; EXHAUST_SIGMA=0.12
MODES=['Still','Uniform','Recirculation','Unsteady']
DIR=np.array([[1,0],[-1,0],[0,1],[0,-1]])
class Simulation:
    def __init__(self):
        self.x=(np.arange(N)+.5)*DX
        self.X,self.Y=np.meshgrid(self.x,self.x)
        self.xf,self.yf=np.meshgrid(np.arange(N+1)*DX,np.arange(N+1)*DX)
        self.mode='Recirculation'; self.gain=.5; self.wind=.15
        self.suction=SUCTION
        k=np.arange(N); eig=-4*np.sin(np.pi*k/(2*N))**2/DX**2
        self.eig=eig[:,None]+eig[None,:]; self.eig[0,0]=1
        self._flow_key=None
        self.reset()
    def reset(self):
        self.C=np.zeros((N,N)); self.t=0.; self.source=np.array([1.2,3.6]); self.robot=np.array([3.6,1.4]); self.path=[self.robot.copy()]; self.history=[]
        self.velocity=np.zeros(2); self.readings=np.zeros(4); self.reached=False
    def ambient_flow(self):
        if self.mode=='Still': return np.zeros((N,N+1)),np.zeros((N+1,N))
        if self.mode=='Uniform': return np.full((N,N+1),self.wind),np.zeros((N+1,N))
        # Discrete curl of streamfunction: exactly divergence-free on cell faces.
        p=self.wind*L/np.pi*np.sin(np.pi*self.xf/L)*np.sin(np.pi*self.yf/L)
        if self.mode=='Unsteady':
            for m,n,a,w in [(2,1,.22,.13),(1,2,.19,.19),(3,2,.10,.31)]:
                p+=self.wind*L/np.pi*a*np.sin(w*self.t+m)*np.sin(m*np.pi*self.xf/L)*np.sin(n*np.pi*self.yf/L)
        p[[0,-1],:]=0; p[:,[0,-1]]=0
        return np.diff(p,axis=0)/DX,-np.diff(p,axis=1)/DX
    def aspiration_flow(self):
        # Depth-averaged volumetric flux per unit room height (m^2/s).
        # Four equal intake sinks, balanced by a diffuse central exhaust.
        key=(*self.robot,self.suction)
        if key==self._flow_key: return self._asp_flow
        def kernel(center,sigma):
            g=np.exp(-((self.X-center[0])**2+(self.Y-center[1])**2)/(2*sigma**2))
            return g/(g.sum()*DX**2)
        intake=sum(kernel(self.robot+INLET*d,PORT_SIGMA) for d in DIR)*self.suction
        exhaust=4*self.suction*kernel(self.robot,EXHAUST_SIGMA)
        rhs=exhaust-intake
        # Neumann Poisson solve: div(grad(phi)) = exhaust - intake.
        coeff=dctn(rhs,type=2,norm='ortho')/self.eig; coeff[0,0]=0
        phi=idctn(coeff,type=2,norm='ortho')
        u=np.zeros((N,N+1)); v=np.zeros((N+1,N))
        u[:,1:-1]=np.diff(phi,axis=1)/DX; v[1:-1,:]=np.diff(phi,axis=0)/DX
        self._flow_key=key; self._asp_flow=(u,v,intake,exhaust)
        return self._asp_flow
    def flow(self):
        u,v=self.ambient_flow(); a,b,_,_=self.aspiration_flow()
        return u+a,v+b
    def sample(self,points):
        q=np.clip(points/DX-.5,0,N-1); i=np.minimum(q[:,0].astype(int),N-2); j=np.minimum(q[:,1].astype(int),N-2)
        a=q[:,0]-i; b=q[:,1]-j; c=self.C
        return (1-a)*(1-b)*c[j,i]+a*(1-b)*c[j,i+1]+(1-a)*b*c[j+1,i]+a*b*c[j+1,i+1]
    def sense(self):
        # Equal flow channels, weighted footprint in front of each inlet.
        distances=np.linspace(INLET,INLET+REACH,6); lateral=np.linspace(-.07,.07,5)
        r,s=np.meshgrid(distances,lateral); weights=np.exp(-(r-INLET)/.15-(s/.06)**2).ravel(); weights/=weights.sum()
        return np.array([self.sample(self.robot+r.ravel()[:,None]*d+s.ravel()[:,None]*np.array([-d[1],d[0]]))@weights for d in DIR])
    def step(self,move=True):
        if self.reached and move: return
        u,v=self.flow(); c=self.C
        _,_,intake,exhaust=self.aspiration_flow()
        # Aspirated odor is mixed and returned through the exhaust.
        mixed=(intake*c).sum()*DX**2/(4*self.suction) if self.suction else 0.
        # Conservative upwind advection and centered diffusion. Clean inflow,
        # zero diffusive wall flux; uniform mode has open left/right boundaries.
        l=np.pad(c,((0,0),(1,0)),mode='edge'); r=np.pad(c,((0,0),(0,1)),mode='edge')
        b=np.pad(c,((1,0),(0,0)),mode='edge'); a=np.pad(c,((0,1),(0,0)),mode='edge')
        fx=u*np.where(u>=0,l,r)-D*(r-l)/DX; fy=v*np.where(v>=0,b,a)-D*(a-b)/DX
        if self.mode=='Uniform': fx[:,0]=0 # clean incoming air
        rate=(np.max(abs(u))+np.max(abs(v)))/DX+4*D/DX**2+DECAY+np.max(intake)
        if DT*rate>1: raise RuntimeError('Unstable timestep: reduce DT')
        src=np.exp(-((self.X-self.source[0])**2+(self.Y-self.source[1])**2)/(2*.05**2)); src*=EMISSION/(src.sum()*DX**2)
        self.C=c+DT*(-(np.diff(fx,axis=1)+np.diff(fy,axis=0))/DX+src-DECAY*c-intake*c+exhaust*mixed)
        self.t+=DT; self.readings=self.sense(); p,m,q,n=self.readings
        # Raw proportional gradient ascent, no knowledge of source position.
        self.velocity=self.gain*np.array([p-m,q-n]) if move else np.zeros(2)
        speed=np.linalg.norm(self.velocity)
        if speed>MAX_SPEED: self.velocity*=MAX_SPEED/speed
        old=self.robot.copy(); self.robot=np.clip(self.robot+DT*self.velocity,INLET+REACH+.08,L-INLET-REACH-.08); self.velocity=(self.robot-old)/DT
        self.readings=self.sense()
        if move and np.linalg.norm(self.robot-self.source)<=CAPTURE_RADIUS:
            self.reached=True; self.velocity[:]=0
    def record(self):
        self.path.append(self.robot.copy()); self.history.append([self.t,*self.readings])

def run():
    s=Simulation(); fig=plt.figure(figsize=(13,8)); gs=fig.add_gridspec(1,2,width_ratios=[2.1,1],left=.065,right=.96,bottom=.25,top=.94,wspace=.32)
    room=fig.add_subplot(gs[0]); signals=fig.add_subplot(gs[1]); im=room.imshow(s.C,origin='lower',extent=[0,L,0,L],cmap='magma',vmin=0,vmax=1)
    fig.colorbar(im,ax=room,fraction=.035,pad=.025,label='Concentration (a.u.)')
    room.set(xlabel='x (m)',ylabel='y (m)',title='Room — left-click source; right-click robot'); room.set_aspect('equal')
    source,=room.plot(*s.source,'w*',ms=13,label='Source',zorder=7)
    robot=Circle(s.robot,BODY_RADIUS,facecolor='#164d58',edgecolor='cyan',lw=1.5,label='Robot',zorder=8); room.add_patch(robot)
    path,=room.plot([],[],'c-',lw=1,alpha=.7,zorder=5)
    ports=[Circle(s.robot+INLET*d,.025,facecolor='white',edgecolor='#164d58',zorder=9) for d in DIR]
    for port in ports: room.add_patch(port)
    inlet_arrows=[FancyArrowPatch(s.robot+.25*d,s.robot+.16*d,arrowstyle='-|>',mutation_scale=8,color='cyan',alpha=.8,zorder=9) for d in DIR]
    for arrow in inlet_arrows: room.add_patch(arrow)
    xx,yy=np.meshgrid(s.x[::7],s.x[::7]); arrows=room.quiver(xx,yy,np.zeros_like(xx),np.zeros_like(yy),color='white',alpha=.32,angles='xy',scale_units='xy',scale=1.2,width=.0025,pivot='mid')
    # A sparse local refinement of the same TOTAL velocity field, not a second flow.
    local=np.linspace(-.32,.32,5); lx,ly=np.meshgrid(local,local)
    offsets=np.column_stack([lx.ravel(),ly.ravel()]); offsets=offsets[np.linalg.norm(offsets,axis=1)>.19]
    near=room.quiver(offsets[:,0],offsets[:,1],np.zeros(len(offsets)),np.zeros(len(offsets)),color='cyan',angles='xy',scale_units='xy',scale=2.5,width=.0018,alpha=.32,pivot='mid',zorder=4)
    room.set_xlim(0,L); room.set_ylim(0,L); room.set_autoscale_on(False)
    room.legend(loc='upper right',fontsize=8)
    lines=[signals.plot([],[],label=k)[0] for k in ['+x','-x','+y','-y']]; signals.set(title='Four aspiration signals',ylabel='Concentration (a.u.)',xlabel='Time (s)'); signals.legend(ncol=2); signals.grid(alpha=.2)
    status=fig.text(.5,.015,'',ha='center'); paused=False
    radio=RadioButtons(fig.add_axes([.025,.04,.15,.17]),MODES,active=2)
    gain=Slider(fig.add_axes([.29,.18,.45,.022]),'Gain',0,5,valinit=s.gain)
    wind=Slider(fig.add_axes([.29,.14,.45,.022]),'Wind (m/s)',0,.25,valinit=s.wind)
    suction=Slider(fig.add_axes([.29,.10,.45,.022]),'Suction (m²/s)',0,.02,valinit=s.suction)
    speed=Slider(fig.add_axes([.29,.06,.45,.022]),'Speed (steps/frame)',1,100,valinit=12,valstep=1)
    pause=Button(fig.add_axes([.81,.16,.15,.04]),'Pause'); reset=Button(fig.add_axes([.81,.105,.15,.04]),'Reset'); warm=Button(fig.add_axes([.81,.05,.15,.04]),'Pre-fill 60 s')
    radio.on_clicked(lambda val:setattr(s,'mode',val)); gain.on_changed(lambda val:setattr(s,'gain',val)); wind.on_changed(lambda val:setattr(s,'wind',val)); suction.on_changed(lambda val:setattr(s,'suction',val))
    def toggle(_):
        nonlocal paused
        if s.reached: return
        paused=not paused; pause.label.set_text('Resume' if paused else 'Pause')
    def restart(_):
        nonlocal paused
        s.reset(); paused=False; pause.label.set_text('Pause')
    pause.on_clicked(toggle); reset.on_clicked(restart)
    def prefill(_):
        if s.reached: return
        for i in range(round(60/DT)):
            s.step(move=False)
            if i%12==0: s.record()
        s.record()
    warm.on_clicked(prefill)
    def click(e):
        nonlocal paused
        if e.inaxes!=room or e.xdata is None: return
        if e.button==1: s.source=np.clip([e.xdata,e.ydata],.48,L-.48)
        elif e.button==3: s.robot=np.clip([e.xdata,e.ydata],.48,L-.48); s.path=[s.robot.copy()]
        else: return
        if s.reached: s.reached=False; paused=True; pause.label.set_text('Resume')
    fig.canvas.mpl_connect('button_press_event',click)
    def update(_):
        if not paused and not s.reached:
            for _ in range(int(speed.val)):
                s.step()
                if s.reached: break
            s.record()
        im.set_data(s.C); im.set_clim(0,max(.01,float(s.C.max())))
        source.set_data([s.source[0]],[s.source[1]]); robot.center=tuple(s.robot)
        p=np.asarray(s.path); path.set_data(p[:,0],p[:,1]); 
        for d,port,arrow in zip(DIR,ports,inlet_arrows):
            port.center=tuple(s.robot+INLET*d); arrow.set_positions(s.robot+.25*d,s.robot+.16*d)
        u,v=s.flow(); uc=(u[:,:-1]+u[:,1:])/2; vc=(v[:-1,:]+v[1:,:])/2
        mask=(xx-s.robot[0])**2+(yy-s.robot[1])**2<.42**2
        arrows.set_UVC(np.ma.array(uc[::7,::7],mask=mask),np.ma.array(vc[::7,::7],mask=mask)); pts=s.robot+offsets; near.set_offsets(pts)
        ij=np.clip((pts/DX).astype(int),0,N-1); near.set_UVC(uc[ij[:,1],ij[:,0]],vc[ij[:,1],ij[:,0]])
        h=np.asarray(s.history)
        for k,line in enumerate(lines): line.set_data(h[:,0],h[:,k+1]) if len(h) else line.set_data([],[])
        signals.set_xlim(0,max(10,s.t)); signals.relim(); signals.autoscale_view(scalex=False); signals.set_ylim(bottom=0)
        status.set_text(f't = {s.t:.1f} s | '+('SOURCE REACHED — simulation stopped' if s.reached else 'Paused' if paused else 'Running')+f' | {int(speed.val)*DT:.2f} simulated s/frame')
        if s.reached: pause.label.set_text('Reached')
    anim=FuncAnimation(fig,update,interval=50,cache_frame_data=False); plt.show()

def test():
    for mode in MODES:
        s=Simulation(); s.mode=mode; s.wind=.25
        for _ in range(400): s.step()
        assert np.isfinite(s.C).all() and s.C.min()>=-1e-12
        mass=s.C.sum()*DX**2
        assert 0<mass<=EMISSION*s.t+1e-8
        print(mode, 'OK', 'mass=',round(mass,4))
    s=Simulation(); s.mode='Still'; s.C=s.X.copy(); z=s.sense(); assert z[0]>z[1] and abs(z[2]-z[3])<1e-10
    print('Directional sensing OK')
    s=Simulation(); s.mode='Still'; u,v,i,e=s.aspiration_flow()
    assert np.max(abs((np.diff(u,axis=1)+np.diff(v,axis=0))/DX-(e-i)))<1e-9
    s.C[:]=2; s.step(move=False); assert np.isclose(s.C.sum()*DX**2,50+DT*(EMISSION-DECAY*50))
    s.robot=s.source.copy(); s.step(); assert s.reached; t=s.t; s.step(); assert s.t==t
    print('Aspiration continuity, scalar mass balance and automatic stop OK')
if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--test',action='store_true'); args=parser.parse_args(); test() if args.test else run()