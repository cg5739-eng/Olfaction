/* Browser port of odor_sim.py. No external dependencies. */
(function(root){
"use strict";
const L=5,N=100,DX=L/N,DT=.025,D=.005,DECAY=.002,EMISSION=1,INLET=.15,REACH=.25,MAX_SPEED=.2,CAPTURE=.2;
const MODES=["Still","Uniform","Recirculation","Unsteady"], DIR=[[1,0],[-1,0],[0,1],[0,-1]];
const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const field=()=>new Float64Array(N*N);
const basis=field(), inverse=field(), eig=field();
for(let k=0;k<N;k++)for(let j=0;j<N;j++){
  basis[k*N+j]=Math.sqrt((k?2:1)/N)*Math.cos(Math.PI*k*(j+.5)/N);
  inverse[j*N+k]=basis[k*N+j];
  eig[k*N+j]=-4*(Math.sin(Math.PI*k/(2*N))**2+Math.sin(Math.PI*j/(2*N))**2)/DX**2;
}
// Separable orthonormal DCT-II and its inverse; same Neumann eigenvalues as SciPy.
function transform(a,m,tmp,out){
  for(let y=0;y<N;y++)for(let k=0;k<N;k++){
    let s=0, r=y*N, q=k*N;
    for(let j=0;j<N;j++)s+=a[r+j]*m[q+j];
    tmp[k*N+y]=s;
  }
  for(let x=0;x<N;x++)for(let k=0;k<N;k++){
    let s=0,r=x*N,q=k*N;
    for(let j=0;j<N;j++)s+=tmp[r+j]*m[q+j];
    out[k*N+x]=s;
  }
}
function kernel(x,y,sigma){
  const g=field(), gx=new Float64Array(N),gy=new Float64Array(N);
  let sx=0,sy=0;
  for(let i=0;i<N;i++){gx[i]=Math.exp(-(((i+.5)*DX-x)**2)/(2*sigma*sigma));gy[i]=Math.exp(-(((i+.5)*DX-y)**2)/(2*sigma*sigma));sx+=gx[i];sy+=gy[i];}
  const f=1/(sx*sy*DX*DX);
  for(let j=0;j<N;j++)for(let i=0;i<N;i++)g[j*N+i]=gx[i]*gy[j]*f;
  return g;
}
const samples=[];let ws=0;
for(let j=0;j<5;j++)for(let i=0;i<6;i++){
  const r=INLET+REACH*i/5,s=-.07+.14*j/4,w=Math.exp(-(r-INLET)/.15-(s/.06)**2);
  samples.push([r,s,w]);ws+=w;
}
for(const a of samples)a[2]/=ws;
class Simulation {
  constructor(){this.mode="Recirculation";this.gain=.5;this.wind=.15;this.suction=.012;
    this.tmp=field();this.coeff=field();this.phi=field();this.rhs=field();this.next=field();
    this.u=new Float64Array(N*(N+1));this.v=new Float64Array(N*(N+1));
    this.fx=new Float64Array(this.u.length);this.fy=new Float64Array(this.v.length);
    this.reset();}
  reset(){this.C=field();this.t=0;this.source=[1.2,3.6];this.robot=[3.6,1.4];this.reached=false;this.velocity=[0,0];this.readings=[0,0,0,0];this.aspKey=null;this.ambKey=null;this.srcKey=null;}
  ambient(){
    const key=[this.mode,this.wind,this.mode==="Unsteady"?this.t:0].join(",");
    if(key===this.ambKey)return this.amb;
    const u=new Float64Array(this.u.length),v=new Float64Array(this.v.length);
    if(this.mode==="Uniform")u.fill(this.wind);
    else if(this.mode!=="Still"){
      const p=new Float64Array((N+1)**2),a=this.wind*L/Math.PI;
      const modes=this.mode==="Unsteady"?[[1,1,1],[2,1,.22*Math.sin(.13*this.t+2)],[1,2,.19*Math.sin(.19*this.t+1)],[3,2,.10*Math.sin(.31*this.t+3)]]:[[1,1,1]];
      for(const [m,n,b] of modes){
        const sx=Array.from({length:N+1},(_,i)=>Math.sin(m*Math.PI*i/N));
        const sy=Array.from({length:N+1},(_,j)=>Math.sin(n*Math.PI*j/N));
        for(let j=1;j<N;j++)for(let i=1;i<N;i++)p[j*(N+1)+i]+=a*b*sx[i]*sy[j];
      }
      for(let j=0;j<N;j++)for(let i=0;i<=N;i++)u[j*(N+1)+i]=(p[(j+1)*(N+1)+i]-p[j*(N+1)+i])/DX;
      for(let j=0;j<=N;j++)for(let i=0;i<N;i++)v[j*N+i]=-(p[j*(N+1)+i+1]-p[j*(N+1)+i])/DX;
    }
    this.ambKey=key;return this.amb={u,v};
  }
  aspiration(){
    const key=[...this.robot,this.suction].join(",");
    if(key===this.aspKey)return this.asp;
    const u=new Float64Array(this.u.length),v=new Float64Array(this.v.length),intake=field(),exhaust=field(),[x,y]=this.robot,q=this.suction;
    if(q>0){
      for(const [dx,dy] of DIR){const g=kernel(x+INLET*dx,y+INLET*dy,.055);for(let k=0;k<g.length;k++)intake[k]+=q*g[k];}
      const g=kernel(x,y,.12);
      for(let k=0;k<g.length;k++){exhaust[k]=4*q*g[k];this.rhs[k]=exhaust[k]-intake[k];}
      transform(this.rhs,basis,this.tmp,this.coeff);
      this.coeff[0]=0;
      for(let k=1;k<this.coeff.length;k++)this.coeff[k]/=eig[k];
      transform(this.coeff,inverse,this.tmp,this.phi);
      for(let j=0;j<N;j++)for(let i=1;i<N;i++)u[j*(N+1)+i]=(this.phi[j*N+i]-this.phi[j*N+i-1])/DX;
      for(let j=1;j<N;j++)for(let i=0;i<N;i++)v[j*N+i]=(this.phi[j*N+i]-this.phi[(j-1)*N+i])/DX;
    }
    this.aspKey=key;return this.asp={u,v,intake,exhaust};
  }
  flow(){
    const a=this.ambient(),b=this.aspiration();
    for(let k=0;k<this.u.length;k++){this.u[k]=a.u[k]+b.u[k];this.v[k]=a.v[k]+b.v[k];}
    return {u:this.u,v:this.v,intake:b.intake,exhaust:b.exhaust};
  }
  sample(x,y){
    const q=clamp(x/DX-.5,0,N-1),p=clamp(y/DX-.5,0,N-1);
    const i=Math.min(Math.floor(q),N-2),j=Math.min(Math.floor(p),N-2),a=q-i,b=p-j,c=this.C,k=j*N+i;
    return (1-a)*(1-b)*c[k]+a*(1-b)*c[k+1]+(1-a)*b*c[k+N]+a*b*c[k+N+1];
  }
  sense(){const [x,y]=this.robot;return DIR.map(([dx,dy])=>samples.reduce((sum,[r,s,w])=>sum+w*this.sample(x+r*dx-s*dy,y+r*dy+s*dx),0));}
  emission(){
    const key=this.source.join(",");
    if(key!==this.srcKey){this.src=kernel(...this.source,.05);this.srcKey=key;}
    return this.src;
  }
  step(move=true){
    if(this.reached&&move)return;
    const {u,v,intake,exhaust}=this.flow(),c=this.C,src=this.emission(),fx=this.fx,fy=this.fy;
    let mixed=0,um=0,vm=0,im=0;
    for(let k=0;k<c.length;k++){mixed+=intake[k]*c[k];im=Math.max(im,intake[k]);}
    mixed=this.suction?mixed*DX*DX/(4*this.suction):0;
    for(let k=0;k<u.length;k++){um=Math.max(um,Math.abs(u[k]));vm=Math.max(vm,Math.abs(v[k]));}
    if(DT*((um+vm)/DX+4*D/DX**2+DECAY+im)>1)throw Error("Unstable numerical timestep.");
    for(let j=0;j<N;j++)for(let i=0;i<=N;i++){
      const k=j*(N+1)+i,l=c[j*N+Math.max(0,i-1)],r=c[j*N+Math.min(N-1,i)];
      fx[k]=(i===0&&this.mode==="Uniform")?0:u[k]*(u[k]>=0?l:r)-D*(r-l)/DX;
    }
    for(let j=0;j<=N;j++)for(let i=0;i<N;i++){
      const k=j*N+i,b=c[Math.max(0,j-1)*N+i],a=c[Math.min(N-1,j)*N+i];
      fy[k]=v[k]*(v[k]>=0?b:a)-D*(a-b)/DX;
    }
    for(let j=0;j<N;j++)for(let i=0;i<N;i++){
      const k=j*N+i,f=j*(N+1)+i;
      this.next[k]=c[k]+DT*(-(fx[f+1]-fx[f]+fy[k+N]-fy[k])/DX+EMISSION*src[k]-DECAY*c[k]-intake[k]*c[k]+exhaust[k]*mixed);
    }
    this.C=this.next;this.next=c;this.t+=DT;
    const [p,m,q,n]=this.sense();
    let vx=move?this.gain*(p-m):0,vy=move?this.gain*(q-n):0,sp=Math.hypot(vx,vy);
    if(sp>MAX_SPEED){vx*=MAX_SPEED/sp;vy*=MAX_SPEED/sp;}
    const [x,y]=this.robot;
    this.robot=[clamp(x+DT*vx,.48,L-.48),clamp(y+DT*vy,.48,L-.48)];
    this.velocity=[(this.robot[0]-x)/DT,(this.robot[1]-y)/DT];this.readings=this.sense();
    if(move&&Math.hypot(this.robot[0]-this.source[0],this.robot[1]-this.source[1])<=CAPTURE){this.reached=true;this.velocity=[0,0];}
  }
}
root.Olfaction={Simulation,L,N,DX,DT,DIR,MODES,INLET,CAPTURE};
})(typeof globalThis!=="undefined"?globalThis:this);
