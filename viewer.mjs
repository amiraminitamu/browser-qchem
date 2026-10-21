// SPDX-License-Identifier: AGPL-3.0-only
import {ANGSTROM_PER_BOHR} from './molecule.mjs?v=e933e699b7cf09943ce2';
import {inferBonds} from './scan.mjs?v=e933e699b7cf09943ce2';
import {COORDINATES,coordinateValue,add,sub,scale,cross,normalize,rotateVector} from './coordinates.mjs?v=e933e699b7cf09943ce2';
const vector=a=>({x:a[0],y:a[1],z:a[2]});
const AMBER='#d96b27';
export class MoleculeViewer {
  constructor(element,onSelect,onBondSelect=()=>{}) {
    this.element=element;this.onSelect=onSelect;this.onBondSelect=onBondSelect;this.model=null;this.labels=true;
    this.pair=[];this.moving=[];this.picking='none';this.coordinate='bond';this.shapes=[];this.ownLabels=[];
    try {
      if(!globalThis.$3Dmol)throw new Error('Viewer library could not be loaded.');
      this.viewer=globalThis.$3Dmol.createViewer(element,{backgroundColor:'#f0ece3',antialias:true,minimumZoomToDistance:2.6});
      if(!this.viewer)throw new Error('WebGL is unavailable.');
      this.viewer.setProjection('orthographic');
      new ResizeObserver(()=>{this.viewer.resize();if(this.model)this.fit();}).observe(element);
    }catch(error){element.textContent=`3D viewer unavailable: ${error.message} Calculations and atom-index controls still work.`;element.classList.add('viewer-fallback');}
  }
  setMolecule(molecule,{reset=false}={}) {
    this.molecule=molecule;if(!this.viewer)return;
    const oldView=this.model&&!reset?this.viewer.getView():null;
    this.clearOverlays();this.viewer.removeAllModels();
    this.bonds=inferBonds(molecule);
    const atoms=molecule.atoms.map((atom,index)=>({elem:atom.symbol,index,serial:index+1,...vector(scale(atom.xyz,ANGSTROM_PER_BOHR)),bonds:[],bondOrder:[]}));
    for(const [a,b] of this.bonds){atoms[a].bonds.push(b);atoms[a].bondOrder.push(1);atoms[b].bonds.push(a);atoms[b].bondOrder.push(1);}
    this.model=this.viewer.addModel();this.model.addAtoms(atoms);this.bindPicking();this.style();
    if(oldView)this.viewer.setView(oldView);else{if(!this.oriented){this.viewer.rotate(65,'y');this.viewer.rotate(12,'x');this.oriented=true;}this.fit();}
    this.viewer.render();this.element.dataset.atomCount=String(atoms.length);
  }
  setPickingMode(mode,coordinate='bond') {
    const changed=this.picking!==mode||this.coordinate!==coordinate;
    this.picking=mode;this.coordinate=coordinate;this.element.dataset.picking=mode;
    if(changed&&this.model){this.bindPicking();this.style();}
  }
  bindPicking() {
    if(!this.viewer||!this.model)return;
    const active=this.picking==='internal-coordinate';
    this.viewer.setClickable({},active,atom=>{if(this.picking==='internal-coordinate')this.onSelect(atom.index);});
    this.viewer.setHoverable({},active,atom=>this.hoverAtom(atom),()=>this.clearHover());
  }
  clearHover() {
    if(!this.viewer)return;
    if(this.hoverShape)this.viewer.removeShape(this.hoverShape);
    if(this.hoverLabel)this.viewer.removeLabel(this.hoverLabel);
    this.hoverShape=null;this.hoverLabel=null;this.viewer.render();
  }
  hoverAtom(atom) {
    if(this.picking!=='internal-coordinate')return;
    this.clearHover();
    this.hoverShape=this.viewer.addSphere({center:atom,radius:.5,color:AMBER,wireframe:true,opacity:.5});
    const next=Math.min(this.pair.length+1,COORDINATES[this.coordinate].arity);
    this.hoverLabel=this.viewer.addLabel(`${atom.elem}${atom.index+1} · Click to assign Atom ${this.pair.length>=COORDINATES[this.coordinate].arity?1:next}`,{position:atom,fontSize:11,fontColor:'#282e38',backgroundColor:'#fff',backgroundOpacity:.96,borderRadius:4,inFront:true,screenOffset:{x:15,y:-22}});
    this.viewer.render();
  }
  clearOverlays() {
    if(!this.viewer)return;
    this.clearHover();
    for(const shape of this.shapes)this.viewer.removeShape(shape);
    for(const label of this.ownLabels)this.viewer.removeLabel(label);
    this.shapes=[];this.ownLabels=[];
  }
  shape(method,spec){const result=this.viewer[method](spec);this.shapes.push(result);return result;}
  label(text,position,selected=false) {
    this.ownLabels.push(this.viewer.addLabel(text,{position:vector(position),fontSize:selected?12:11,fontColor:selected?'#fff':'#515b68',backgroundColor:selected?'#ad4c16':'#ffffff',backgroundOpacity:selected?1:.9,borderRadius:4,inFront:true,screenOffset:{x:14,y:14}}));
  }
  style() {
    if(!this.viewer||!this.model)return;
    this.clearOverlays();
    this.viewer.setStyle({},{stick:{radius:.12,colorscheme:'Jmol'},sphere:{scale:.32,colorscheme:'Jmol'}});
    const active=this.picking==='internal-coordinate',positions=this.molecule.atoms.map(atom=>scale(atom.xyz,ANGSTROM_PER_BOHR));
    if(this.labels)this.molecule.atoms.forEach((atom,i)=>{if(!active||!this.pair.includes(i))this.label(`${atom.symbol}${i+1}`,positions[i]);});
    if(active) {
      this.pair.forEach((index,order)=>{if(!positions[index])return;this.shape('addSphere',{center:vector(positions[index]),radius:.47,color:AMBER,wireframe:true,opacity:.7});this.label(`${order+1}: ${this.molecule.atoms[index].symbol}${index+1}`,positions[index],true);});
      if(this.coordinate==='bond')for(const [a,b] of this.bonds) {
        const start=add(positions[a],scale(sub(positions[b],positions[a]),.3)),end=add(positions[a],scale(sub(positions[b],positions[a]),.7));
        const shape=this.shape('addCylinder',{start:vector(start),end:vector(end),radius:.16,color:AMBER,opacity:0,clickable:true,hoverable:true,
          callback:()=>{if(this.picking==='internal-coordinate')this.onBondSelect(a,b);},
          hover_callback:()=>{if(this.picking==='internal-coordinate'){shape.updateStyle({opacity:.7});this.viewer.render();}},
          unhover_callback:()=>{shape.updateStyle({opacity:0});this.viewer.render();}});
      }
      this.measurement(positions);
    }
    this.viewer.render();
  }
  measurement(positions) {
    const ids=this.pair,points=ids.map(i=>positions[i]);if(points.some(p=>!p))return;
    for(let i=1;i<points.length;i++)this.shape('addLine',{start:vector(points[i-1]),end:vector(points[i]),color:AMBER,linewidth:2,dashed:true});
    const type=({2:'bond',3:'angle',4:'dihedral'})[points.length];if(!type)return;
    try {
      const measured=coordinateValue(this.molecule,type,ids),text=`${measured.toFixed(type==='bond'?4:2).replace('-', '−')} ${COORDINATES[type].unit}`;
      if(type==='bond'){this.label(text,scale(add(points[0],points[1]),.5));return;}
      if(type==='angle') {
        const origin=points[1],u=normalize(sub(points[0],origin)),v=normalize(sub(points[2],origin));let axis=cross(u,v);
        if(Math.hypot(...axis)<1e-8)axis=cross(u,Math.abs(u[0])<.9?[1,0,0]:[0,1,0]);axis=normalize(axis);
        const radius=Math.min(.65,Math.hypot(...sub(points[0],origin))*.55,Math.hypot(...sub(points[2],origin))*.55);
        const arc=Array.from({length:25},(_,i)=>add(origin,scale(rotateVector(u,axis,measured*Math.PI/180*i/24),radius)));
        this.fan(origin,arc);this.label(text,arc[12]);return;
      }
      this.plane(points.slice(0,3));this.plane(points.slice(1));
      const axis=normalize(sub(points[2],points[1])),origin=scale(add(points[1],points[2]),.5),arm=sub(points[0],points[1]);
      const dot=arm.reduce((s,x,i)=>s+x*axis[i],0),u=normalize(sub(arm,scale(axis,dot)));
      const angle=measured*Math.PI/180,arc=Array.from({length:25},(_,i)=>add(origin,scale(rotateVector(u,axis,angle*i/24),.6)));
      for(let i=1;i<arc.length;i++)this.shape('addLine',{start:vector(arc[i-1]),end:vector(arc[i]),color:AMBER,linewidth:2});
      if(Math.abs(angle)>.01)this.shape('addArrow',{start:vector(arc.at(-3)),end:vector(arc.at(-1)),radius:.025,radiusRatio:2.5,mid:.6,color:AMBER});
      this.label(text,arc[12]);
    }catch{/* A partial/undefined coordinate still shows atom order and connecting lines. */}
  }
  plane(points) {this.shape('addCustom',{vertexArr:points.map(vector),faceArr:[0,1,2,2,1,0],color:AMBER,opacity:.15});}
  fan(origin,arc) {
    const vertices=[origin,...arc],faces=[];for(let i=1;i<vertices.length-1;i++)faces.push(0,i,i+1,i+1,i,0);
    this.shape('addCustom',{vertexArr:vertices.map(vector),faceArr:faces,color:AMBER,opacity:.18});
    for(let i=1;i<arc.length;i++)this.shape('addLine',{start:vector(arc[i-1]),end:vector(arc[i]),color:AMBER,linewidth:2});
  }
  select(pair,moving=[]) {this.pair=[...pair];this.moving=[...moving];this.style();}
  setLabels(value) {this.labels=value;this.style();}
  fit() {
    if(!this.viewer||!this.model)return;
    this.viewer.zoomTo();const corners=[];
    for(const atom of this.model.selectedAtoms({}))for(const dx of [-.65,.65])for(const dy of [-.65,.65])for(const dz of [-.65,.65])corners.push({x:atom.x+dx,y:atom.y+dy,z:atom.z+dz});
    const points=this.viewer.modelToScreen(corners);
    if(points.length){const width=Math.max(...points.map(p=>p.x))-Math.min(...points.map(p=>p.x)),height=Math.max(...points.map(p=>p.y))-Math.min(...points.map(p=>p.y));const zoom=Math.min(this.element.clientWidth*.72/width,this.element.clientHeight*.72/height);if(Number.isFinite(zoom)&&zoom>0)this.viewer.zoom(Math.min(8,Math.max(.1,zoom)));}
    this.viewer.render();
  }
  reset(){this.fit();}
  png(){return this.viewer?.pngURI();}
}
