/* Canvas UI for the browser simulation. */
(function(){
"use strict";
const {Simulation,L,N,DX,DT,DIR}=Olfaction;
const s=new Simulation(),$=id=>document.getElementById(id);
const map=$("map"),chart=$("chart"),mctx=map.getContext("2d"),cctx=chart.getContext("2d");
const buffer=document.createElement("canvas");buffer.width=buffer.height=N;
const bctx=buffer.getContext("2d"),pixels=bctx.createImageData(N,N);
const colors=["#ff9e77","#69d8ff","#a4e2a0","#d2a5ff"];
const stops=[[8,6,23],[48,16,75],[118,33,107],[182,59,96],[239,128,81],[252,246,189]];
const palette=Array.from({length:256},(_,i)=>{let t=i/255*(stops.length-1),k=Math.min(stops.length-2,Math.floor(t)),f=t-k;return stops[k].map((v,j)=>Math.round(v+(stops[k+1][j]-v)*f));});
let running=false,warming=0,warmTotal=0,afterWarm=false,failure="",dirty=true,started=false;
let history=[],path=[],recordInterval=.15,lastRecord=-Infinity,geom=null,lastDraw=0;
const clamp=(v,a,b)=>Math.max(a,Math.min(b,v));
const fmt=v=>v===0?"0":Math.abs(v)<.001?v.toExponential(1):v.toPrecision(3);
function record(force=false){
  if(!force&&s.t-lastRecord<recordInterval-1e-9)return;
  if(history.length&&Math.abs(history[history.length-1][0]-s.t)<1e-10){history[history.length-1]=[s.t,...s.readings];}
  else history.push([s.t,...s.readings]);
  lastRecord=s.t;
  if(history.length>16000){history=history.filter((_,i)=>i%2===0||i===history.length-1);recordInterval*=2;}
  const prev=path[path.length-1];if(!prev||Math.hypot(s.robot[0]-prev[0],s.robot[1]-prev[1])>.004||force)path.push([...s.robot]);
  if(path.length>16000)path=path.filter((_,i)=>i%2===0||i===path.length-1);
}
function resize(canvas,ctx){
  const rect=canvas.getBoundingClientRect(),w=Math.max(1,rect.width),h=Math.max(1,rect.height),dpr=Math.min(2,window.devicePixelRatio||1);
  if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}
  ctx.setTransform(dpr,0,0,dpr,0,0);return {w,h};
}
function line(ctx,x1,y1,x2,y2,color,width=1){ctx.strokeStyle=color;ctx.lineWidth=width;ctx.beginPath();ctx.moveTo(x1,y1);ctx.lineTo(x2,y2);ctx.stroke();}
function arrow(ctx,x,y,dx,dy,color,scale,maxlen){
  let ax=dx*scale,ay=-dy*scale,m=Math.hypot(ax,ay);
  if(m<1.1)return;if(m>maxlen){ax*=maxlen/m;ay*=maxlen/m;m=maxlen;}
  const x0=x-ax/2,y0=y-ay/2,x1=x+ax/2,y1=y+ay/2,h=Math.min(4,m*.35),ux=ax/m,uy=ay/m;
  ctx.strokeStyle=color;ctx.lineWidth=.8;ctx.beginPath();ctx.moveTo(x0,y0);ctx.lineTo(x1,y1);
  ctx.moveTo(x1-h*ux+h*.55*uy,y1-h*uy-h*.55*ux);ctx.lineTo(x1,y1);ctx.lineTo(x1-h*ux-h*.55*uy,y1-h*uy+h*.55*ux);ctx.stroke();
}
function circle(ctx,x,y,r,fill,stroke,width=1){ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);if(fill){ctx.fillStyle=fill;ctx.fill();}if(stroke){ctx.strokeStyle=stroke;ctx.lineWidth=width;ctx.stroke();}}
function drawMap(){
  const {w,h}=resize(map,mctx),ctx=mctx;
  const size=Math.max(1,Math.min(w-56,h-48)),left=(w-size)/2+8,top=(h-size)/2-7,px=size/L;
  geom={left,top,size};const X=x=>left+x*px,Y=y=>top+(L-y)*px;
  ctx.clearRect(0,0,w,h);
  let peak=.01;for(const c of s.C)peak=Math.max(peak,c);
  for(let j=0;j<N;j++)for(let i=0;i<N;i++){
    const k=4*((N-1-j)*N+i),rgb=palette[clamp(Math.floor(255*s.C[j*N+i]/peak),0,255)];
    pixels.data[k]=rgb[0];pixels.data[k+1]=rgb[1];pixels.data[k+2]=rgb[2];pixels.data[k+3]=255;
  }
  bctx.putImageData(pixels,0,0);ctx.imageSmoothingEnabled=true;ctx.drawImage(buffer,left,top,size,size);
  ctx.save();ctx.beginPath();ctx.rect(left,top,size,size);ctx.clip();
  // All plotted arrows sample the same total flow used by the transport solver.
  const {u,v}=s.flow();
  function flowAt(x,y){const i=clamp(Math.floor(x/DX),0,N-1),j=clamp(Math.floor(y/DX),0,N-1),k=j*(N+1)+i;return [(u[k]+u[k+1])/2,(v[j*N+i]+v[(j+1)*N+i])/2];}
  for(let j=0;j<N;j+=8)for(let i=0;i<N;i+=8){
    const x=(i+.5)*DX,y=(j+.5)*DX;
    if(Math.hypot(x-s.robot[0],y-s.robot[1])<.42)continue;
    arrow(ctx,X(x),Y(y),...flowAt(x,y),"rgba(240,244,255,.32)",px*.8,px*.3);
  }
  for(let j=-2;j<=2;j++)for(let i=-2;i<=2;i++){
    const ox=.16*i,oy=.16*j;if(Math.hypot(ox,oy)<.19)continue;
    const x=s.robot[0]+ox,y=s.robot[1]+oy;
    arrow(ctx,X(x),Y(y),...flowAt(x,y),"rgba(89,220,228,.30)",px*.45,px*.13);
  }
  if(path.length>1){ctx.strokeStyle="rgba(89,220,228,.8)";ctx.lineWidth=1.4;ctx.beginPath();path.forEach(([x,y],i)=>i?ctx.lineTo(X(x),Y(y)):ctx.moveTo(X(x),Y(y)));ctx.stroke();}
  const [sx,sy]=s.source;ctx.beginPath();
  for(let i=0;i<10;i++){const a=-Math.PI/2+i*Math.PI/5,r=(i%2?.42:1)*Math.max(7,px*.085);const x=X(sx)+r*Math.cos(a),y=Y(sy)+r*Math.sin(a);i?ctx.lineTo(x,y):ctx.moveTo(x,y);}
  ctx.closePath();ctx.fillStyle="#fff";ctx.fill();ctx.strokeStyle="#111827";ctx.lineWidth=.8;ctx.stroke();
  const [rx,ry]=s.robot;
  circle(ctx,X(rx),Y(ry),.15*px,"#164d58","#59dce4",1.5);
  for(const [dx,dy] of DIR){
    circle(ctx,X(rx+.15*dx),Y(ry+.15*dy),Math.max(1.8,.023*px),"#fff","#164d58",.7);
    arrow(ctx,X(rx+.215*dx),Y(ry+.215*dy),-dx,-dy,"rgba(89,220,228,.85)",px*.09,px*.09);
  }
  ctx.restore();ctx.strokeStyle="#667484";ctx.lineWidth=1;ctx.strokeRect(left,top,size,size);
  ctx.font="11px system-ui";ctx.fillStyle="#a8b5c6";ctx.textAlign="center";ctx.textBaseline="top";
  for(let i=0;i<=5;i++){ctx.fillText(String(i),X(i),top+size+6);line(ctx,X(i),top+size,X(i),top+size+3,"#667484");}
  ctx.fillText("x (m)",left+size/2,top+size+21);
  ctx.textAlign="right";ctx.textBaseline="middle";
  for(let i=0;i<=5;i++){ctx.fillText(String(i),left-6,Y(i));line(ctx,left-3,Y(i),left,Y(i),"#667484");}
  ctx.save();ctx.translate(left-27,top+size/2);ctx.rotate(-Math.PI/2);ctx.textAlign="center";ctx.fillText("y (m)",0,0);ctx.restore();
  $("scaleMax").textContent=fmt(peak)+" a.u.";
}
function drawChart(){
  const {w,h}=resize(chart,cctx),ctx=cctx,l=47,r=10,t=12,b=30,pw=Math.max(1,w-l-r),ph=Math.max(1,h-t-b);
  ctx.clearRect(0,0,w,h);
  let ymax=.01;for(const row of history)for(let i=1;i<=4;i++)ymax=Math.max(ymax,row[i]);
  ymax*=1.05;const xmax=Math.max(10,s.t),X=x=>l+pw*x/xmax,Y=y=>t+ph-ph*y/ymax;
  ctx.font="11px system-ui";ctx.fillStyle="#a8b5c6";ctx.textAlign="right";ctx.textBaseline="middle";
  for(let i=0;i<=4;i++){const v=ymax*i/4,y=Y(v);line(ctx,l,y,l+pw,y,"#303d4c");ctx.fillText(fmt(v),l-5,y);}
  ctx.textAlign="center";ctx.textBaseline="top";
  for(let i=0;i<=4;i++){const v=xmax*i/4,x=X(v);line(ctx,x,t,x,t+ph,"#303d4c");ctx.fillText(v.toFixed(xmax<20?1:0),x,t+ph+5);}
  ctx.fillText("Time (s)",l+pw/2,h-13);
  ctx.save();ctx.beginPath();ctx.rect(l,t,pw,ph);ctx.clip();
  for(let k=0;k<4;k++){ctx.strokeStyle=colors[k];ctx.lineWidth=1.6;ctx.beginPath();history.forEach((a,i)=>i?ctx.lineTo(X(a[0]),Y(a[k+1])):ctx.moveTo(X(a[0]),Y(a[k+1])));ctx.stroke();}
  ctx.restore();
}
function render(){
  drawMap();drawChart();
  $("time").textContent="t = "+s.t.toFixed(2)+" s";
  s.readings.forEach((v,i)=>$("r"+i).textContent=fmt(v));
  const status=$("status"),state=failure?"error":s.reached?"reached":warming?"warming":running?"running":"paused";
  const labels={error:failure,reached:"SOURCE REACHED · stopped",warming:"Pre-fill "+Math.floor(100*(1-warming/warmTotal))+"%",running:"Running",paused:"Paused"};
  if(status.textContent!==labels[state])status.textContent=labels[state];status.dataset.state=state;
  $("play").textContent=s.reached?"Reached":running||warming?"Pause":started?"Resume":"Start";
  $("play").disabled=s.reached||!!failure;
  $("warm").disabled=s.reached||warming>0||!!failure;
  dirty=false;
}
function clear(){
  running=false;warming=0;warmTotal=0;afterWarm=false;failure="";started=false;
  s.reset();history=[];path=[];lastRecord=-Infinity;recordInterval=.15;record(true);dirty=true;
}
function beginWarm(seconds,resume){
  running=false;warming=Math.round(seconds/DT);warmTotal=warming;afterWarm=resume;started=true;dirty=true;
}
function tick(now){
  try{
    if((running||warming)&&!s.reached&&!failure){
      const count=warming?100:Number($("speed").value),deadline=performance.now()+18;
      for(let k=0;k<count;k++){
        const pre=warming>0;s.step(!pre);record(s.reached);
        if(pre){warming--;if(warming===0){running=afterWarm;record(true);break;}}
        if(s.reached){running=false;record(true);break;}
        if(performance.now()>deadline)break;
      }
      dirty=true;
    }
    if(dirty&&(now-lastDraw>25||s.reached||!running&&!warming)){render();lastDraw=now;}
  }catch(err){running=false;warming=0;failure=String(err.message||err);$("status").textContent=failure;$("status").dataset.state="error";$("play").disabled=true;$("warm").disabled=true;dirty=false;}
  requestAnimationFrame(tick);
}
$("play").addEventListener("click",()=>{
  if(s.reached||failure)return;
  if(warming){warming=0;running=false;afterWarm=false;record(true);}
  else{running=!running;started=true;}
  dirty=true;
});
$("reset").addEventListener("click",clear);
$("warm").addEventListener("click",()=>{if(!warming&&!s.reached)beginWarm(60,running);});
$("easy").addEventListener("click",()=>{
  clear();$("mode").value="Still";$("gain").value=".5";$("wind").value=".15";$("suction").value=".012";sync();
  s.robot=[2,3.6];history=[];path=[];lastRecord=-Infinity;record(true);beginWarm(6,true);
});
function sync(){
  s.mode=$("mode").value;s.gain=Number($("gain").value);s.wind=Number($("wind").value);s.suction=Number($("suction").value);
  $("gainVal").textContent=s.gain.toFixed(2);$("windVal").textContent=s.wind.toFixed(3);$("suctionVal").textContent=s.suction.toFixed(3);$("speedVal").textContent=$("speed").value+" steps/frame";dirty=true;
}
for(const name of ["mode","gain","wind","suction","speed"])$(name).addEventListener("input",sync);
function relocate(target,x,y){
  if(failure)return;
  s[target]=[clamp(x,.48,L-.48),clamp(y,.48,L-.48)];
  if(target==="robot"){path=[[...s.robot]];s.readings=s.sense();}
  if(s.reached){s.reached=false;running=false;}
  record(true);dirty=true;
}
function place(e,target){
  if(!geom)return;const box=map.getBoundingClientRect(),x=e.clientX-box.left-geom.left,y=e.clientY-box.top-geom.top;
  if(x<0||y<0||x>geom.size||y>geom.size)return;
  relocate(target,x/geom.size*L,L-y/geom.size*L);
}
map.addEventListener("click",e=>place(e,$("target").value));
map.addEventListener("contextmenu",e=>{e.preventDefault();place(e,"robot");});
map.addEventListener("keydown",e=>{
  const d={ArrowLeft:[-.1,0],ArrowRight:[.1,0],ArrowUp:[0,.1],ArrowDown:[0,-.1]}[e.key];
  if(!d)return;e.preventDefault();const a=$("target").value;relocate(a,s[a][0]+d[0],s[a][1]+d[1]);
});
window.addEventListener("resize",()=>{dirty=true;});
if(typeof ResizeObserver!=="undefined"){const ro=new ResizeObserver(()=>{dirty=true;});ro.observe(map);ro.observe(chart);}
clear();sync();render();requestAnimationFrame(tick);
})();