// SPDX-License-Identifier: AGPL-3.0-only
const NS='http://www.w3.org/2000/svg';
function node(tag,attributes={},text) {
  const element=document.createElementNS(NS,tag);
  for(const [name,value] of Object.entries(attributes)) element.setAttribute(name,value);
  if(text!==undefined) element.textContent=text;
  return element;
}
export function plot(container,points,{xLabel='Coordinate',yLabel='Energy / Eh',onSelect=()=>{},selected=-1,empty='Run a calculation to see the energy curve.'}={}) {
  if(!points.length){container.replaceChildren();const p=document.createElement('p');p.className='empty-chart';p.textContent=empty;container.append(p);return;}
  const valid=points.filter(p=>Number.isFinite(p.y));
  if(!valid.length){container.textContent='No numerical energies are available yet.';return;}
  const width=760,height=260,left=86,right=24,top=24,bottom=50;
  let xmin=Math.min(...points.map(p=>p.x)),xmax=Math.max(...points.map(p=>p.x));
  let ymin=Math.min(...valid.map(p=>p.y)),ymax=Math.max(...valid.map(p=>p.y));
  if(xmin===xmax){xmin-=.5;xmax+=.5;}
  const pad=Math.max((ymax-ymin)*.12,Math.max(Math.abs(ymax),1)*1e-7);
  ymin-=pad;ymax+=pad;
  const x=value=>left+(value-xmin)/(xmax-xmin)*(width-left-right);
  const y=value=>height-bottom-(value-ymin)/(ymax-ymin)*(height-top-bottom);
  const svg=node('svg',{viewBox:`0 0 ${width} ${height}`,role:'img','aria-label':`${yLabel} versus ${xLabel}`});
  for(let t=0;t<=4;t++) {
    const yy=top+t*(height-top-bottom)/4,value=ymax-t*(ymax-ymin)/4;
    svg.append(node('line',{x1:left,y1:yy,x2:width-right,y2:yy,stroke:'#e6edf2'}));
    svg.append(node('text',{x:left-12,y:yy+4,'text-anchor':'end',class:'axis-tick'},Math.abs(ymax-ymin)<.001?value.toExponential(2):value.toFixed(3)));
    const xx=left+t*(width-left-right)/4;
    svg.append(node('text',{x:xx,y:height-bottom+24,'text-anchor':'middle',class:'axis-tick'},(xmin+t*(xmax-xmin)/4).toFixed(2)));
  }
  let path='';let connected=false;
  for(const p of points) {
    if(!Number.isFinite(p.y)||p.failed){connected=false;continue;}
    path+=`${connected?'L':'M'}${x(p.x)},${y(p.y)} `;connected=true;
  }
  svg.append(node('path',{d:path,fill:'none',stroke:'#159b90','stroke-width':2.5}));
  for(const [i,p] of points.entries()) if(Number.isFinite(p.y)) {
    const circle=node('circle',{cx:x(p.x),cy:y(p.y),r:i===selected?6.5:4.5,
      fill:p.failed?'#ed9f3e':i===selected?'#122c3e':'#ffffff',stroke:p.failed?'#c58324':'#159b90','stroke-width':2,
      tabindex:'0',role:'button','aria-label':p.label??`Point ${i+1}`,'data-point':i});
    circle.append(node('title',{},p.label??`${p.x}: ${p.y}`));
    circle.addEventListener('click',()=>onSelect(i));
    circle.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();onSelect(i);}});
    svg.append(circle);
  }
  svg.append(node('text',{x:(left+width-right)/2,y:height-6,'text-anchor':'middle',class:'axis-title'},xLabel));
  svg.append(node('text',{x:18,y:(top+height-bottom)/2,transform:`rotate(-90 18 ${(top+height-bottom)/2})`,'text-anchor':'middle',class:'axis-title'},yLabel));
  container.replaceChildren(svg);
}
