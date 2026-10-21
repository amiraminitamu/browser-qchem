// SPDX-License-Identifier: AGPL-3.0-only
// One immutable figure model drives the interactive plot and all file formats.
import {COORDINATES,coordinateLabel} from './coordinates.mjs';
import {basisLabel} from './basis.mjs';
import {mathText} from './format.mjs';
const xml=value=>String(value).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
export const ENERGY_UNITS=Object.freeze({eh:{factor:1,label:'Eₕ'},kcal:{factor:627.5094740631,label:'kcal·mol⁻¹'},ev:{factor:27.211386245988,label:'eV'},kj:{factor:2625.4996394799,label:'kJ·mol⁻¹'}});
export function pesFigureModel(scan,{unit='kcal',relative=true}={}) {
  const u=ENERGY_UNITS[unit];
  if(!u)throw new Error('Unknown energy unit.');
  const finite=scan.points.filter(p=>Number.isFinite(p.energy));if(!finite.length)throw new Error('No finite scan energies to export.');
  const converged=finite.filter(p=>p.converged),eligible=converged.length?converged:finite;
  const minimum=eligible.reduce((best,p)=>p.energy<best.energy?p:best);
  const source=scan.source??scan.points[0].molecule,spec=scan.spec,meta=COORDINATES[spec.coordinate??'bond'];
  const coordinate=coordinateLabel(source,spec),symbol=spec.coordinate==='bond'||!spec.coordinate?(spec.unit==='bohr'?'bohr':'Å'):'°';
  return {title:`${scan.moleculeName??'Molecular'} potential energy scan`,xLabel:coordinate,yLabel:`${relative?'Relative energy ΔE':'Total energy'} [${u.label}]`,relative,unit:u.label,minimum,
    points:scan.points.map(p=>({x:p.coordinate,y:Number.isFinite(p.energy)?(p.energy-(relative?minimum.energy:0))*u.factor:null,failed:!p.converged,label:mathText(`Point ${p.index+1}: ${p.coordinate.toFixed(4)} ${symbol}; ${p.energy??'no energy'} Eₕ; ${p.status}`)})),
    metadata:[
      ['Molecule',`${scan.moleculeName??'Custom molecule'}${scan.formula?' · '+scan.formula:''} · ${source.atoms.length} atoms`],
      ['Calculation',`${(spec.method??'rhf').toUpperCase()} / ${basisLabel(spec.basis)}${spec.auxiliaryBasis?' · auxiliary '+spec.auxiliaryBasis:''} · charge ${spec.charge??0} · multiplicity ${spec.multiplicity??1}`],
      ['SCF integrals',spec.scfIntegrals==='df'?`Density fitted (DF-JK) / ${spec.jkfitBasis??'def2-universal-JKFIT'}`:'Exact four-center'],
      ['Coordinate',`${coordinate} · moving atoms ${(spec.moving??[]).map(i=>i+1).join(', ')}`],
      ['Sampled minimum',`${meta.symbol} = ${minimum.coordinate.toFixed(6)} ${symbol}; E = ${minimum.energy.toFixed(12)} Eₕ; trajectory point ${minimum.index+1}${!converged.length?' (unconverged)':''}`],
      ['Scan range',`${spec.start} → ${spec.end} ${symbol} · ${scan.points.length}/${spec.count} points · unrelaxed${scan.cancelled?' · stopped':''}`],
      ['Reference',relative?`Zero at lowest ${converged.length?'converged':'available, unconverged'} sampled energy. Failed/unconverged points are not connected.`:'Absolute electronic total energies. Failed/unconverged points are not connected.']
    ]};
}
export function figureSVG(model) {
  const width=900,chartHeight=470,left=112,right=42,top=58,bottom=78;
  const valid=model.points.filter(p=>Number.isFinite(p.y));
  let xmin=Math.min(...model.points.map(p=>p.x)),xmax=Math.max(...model.points.map(p=>p.x)),ymin=Math.min(...valid.map(p=>p.y)),ymax=Math.max(...valid.map(p=>p.y));
  if(xmin===xmax){xmin-=.5;xmax+=.5;}
  const pad=Math.max((ymax-ymin)*.1,Math.max(Math.abs(ymax),1)*1e-7);ymin-=pad;ymax+=pad;
  const x=v=>left+(v-xmin)/(xmax-xmin)*(width-left-right),y=v=>chartHeight-bottom-(v-ymin)/(ymax-ymin)*(chartHeight-top-bottom);
  const elements=[];
  const text=(px,py,value,extra='')=>elements.push(`<text x="${px}" y="${py}" ${extra}>${xml(mathText(value))}</text>`);
  text(left,27,model.title,'font-size="17" font-weight="600"');
  for(let i=0;i<=5;i++){
    const yy=top+i*(chartHeight-top-bottom)/5,v=ymax-i*(ymax-ymin)/5,xx=left+i*(width-left-right)/5;
    elements.push(`<path d="M${left} ${yy}H${width-right}M${xx} ${top}V${chartHeight-bottom}" fill="none" stroke="#e1e3e6" stroke-width=".7"/>`);
    text(left-12,yy+4,Math.abs(ymax-ymin)<.001?v.toFixed(7):v.toFixed(3),'text-anchor="end" font-family="monospace" font-size="12"');
    text(xx,chartHeight-bottom+24,(xmin+i*(xmax-xmin)/5).toFixed(2),'text-anchor="middle" font-family="monospace" font-size="12"');
  }
  elements.push(`<path d="M${left} ${top}V${chartHeight-bottom}H${width-right}" fill="none" stroke="#303640" stroke-width="1"/>`);
  let path='',connected=false;
  for(const p of model.points){if(!Number.isFinite(p.y)||p.failed){connected=false;continue;}path+=`${connected?'L':'M'}${x(p.x)} ${y(p.y)} `;connected=true;}
  elements.push(`<path d="${path}" fill="none" stroke="#ac4f1d" stroke-width="1.8"/>`);
  for(const p of valid){const xx=x(p.x),yy=y(p.y);elements.push(p.failed?`<path d="M${xx-4} ${yy-4}l8 8m-8 0 8-8" stroke="#ad7223" stroke-width="1.5"/>`:`<circle cx="${xx}" cy="${yy}" r="3.6" fill="#fff" stroke="#ac4f1d" stroke-width="1.5"/>`);}
  text((left+width-right)/2,chartHeight-25,model.xLabel,'text-anchor="middle" font-size="14"');
  const cy=(top+chartHeight-bottom)/2;text(25,cy,model.yLabel,`text-anchor="middle" font-size="14" transform="rotate(-90 25 ${cy})"`);
  // Caption rows wrap independently of viewport width and include immutable job metadata.
  let yy=chartHeight+6;
  elements.push(`<path d="M${left} ${yy-15}H${width-right}" stroke="#bfc3ca"/>`);
  for(const [label,value] of model.metadata){
    const words=mathText(value).split(' '),lines=[''];
    for(const word of words){if((lines.at(-1)+' '+word).length>88)lines.push('');lines[lines.length-1]+=(lines.at(-1)?' ':'')+word;}
    text(left,yy,label.toUpperCase(),'font-size="10" font-weight="600" letter-spacing=".5" fill="#646870"');
    for(const line of lines){text(left+145,yy,line,'font-size="11.5"');yy+=17;}yy+=6;
  }
  const height=yy+16;
  return {width,height,svg:`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}"><rect width="100%" height="100%" fill="#ffffff"/><g fill="#282e38" font-family="Arial,Helvetica,sans-serif">${elements.join('')}</g></svg>`};
}
// A self-contained high-resolution PDF image document, without a remote service.
// Content is generated from scientific data, never from a screenshot of the UI.
export function imagePDF(jpeg,width,height,pageWidth=720) {
  const enc=new TextEncoder(),chunks=[],offsets=[0];let size=0;
  const push=data=>{const bytes=typeof data==='string'?enc.encode(data):data;chunks.push(bytes);size+=bytes.length;};
  const object=(id,body)=>{offsets[id]=size;push(`${id} 0 obj\n${body}\nendobj\n`);};
  const pageHeight=pageWidth*height/width;
  push('%PDF-1.4\n%Lorentzian\n');
  object(1,'<< /Type /Catalog /Pages 2 0 R >>');object(2,'<< /Type /Pages /Kids [3 0 R] /Count 1 >>');
  object(3,`<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${pageWidth} ${pageHeight}] /Resources << /XObject << /Figure 4 0 R >> >> /Contents 5 0 R >>`);
  offsets[4]=size;push(`4 0 obj\n<< /Type /XObject /Subtype /Image /Width ${width} /Height ${height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpeg.length} >>\nstream\n`);push(jpeg);push('\nendstream\nendobj\n');
  const drawing=`q ${pageWidth} 0 0 ${pageHeight} 0 0 cm /Figure Do Q\n`;
  object(5,`<< /Length ${enc.encode(drawing).length} >>\nstream\n${drawing}endstream`);
  object(6,'<< /Title (Lorentzian potential energy scan) /Creator (Lorentzian v0.6) >>');
  const start=size;push('xref\n0 7\n0000000000 65535 f \n');for(let i=1;i<=6;i++)push(`${String(offsets[i]).padStart(10,'0')} 00000 n \n`);
  push(`trailer\n<< /Size 7 /Root 1 0 R /Info 6 0 R >>\nstartxref\n${start}\n%%EOF\n`);
  const result=new Uint8Array(size);let at=0;for(const chunk of chunks){result.set(chunk,at);at+=chunk.length;}return result;
}
export async function exportPESFigure(model,format) {
  if(!['png','jpg','pdf'].includes(format))throw new Error('Unsupported figure format.');
  const figure=figureSVG(model),url=URL.createObjectURL(new Blob([figure.svg],{type:'image/svg+xml'}));
  let blob;
  try{
    const image=new Image();await new Promise((resolve,reject)=>{image.onload=resolve;image.onerror=()=>reject(new Error('Could not render the figure.'));image.src=url;});
    const canvas=document.createElement('canvas');canvas.width=figure.width*3;canvas.height=figure.height*3;
    const context=canvas.getContext('2d');if(!context)throw new Error('Canvas is unavailable.');
    context.fillStyle='#ffffff';context.fillRect(0,0,canvas.width,canvas.height);context.drawImage(image,0,0,canvas.width,canvas.height);
    blob=await new Promise(resolve=>canvas.toBlob(resolve,format==='png'?'image/png':'image/jpeg',format==='pdf'?.98:.95));
    if(!blob)throw new Error('The browser could not encode the figure.');
    if(format==='pdf')blob=new Blob([imagePDF(new Uint8Array(await blob.arrayBuffer()),canvas.width,canvas.height)],{type:'application/pdf'});
  }finally{URL.revokeObjectURL(url);}
  const download=URL.createObjectURL(blob),link=document.createElement('a');link.href=download;link.download=`lorentzian-pes.${format}`;link.click();setTimeout(()=>URL.revokeObjectURL(download),1000);
}
