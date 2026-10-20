// SPDX-License-Identifier: AGPL-3.0-only
import { ANGSTROM_PER_BOHR } from './molecule.mjs';
import { inferBonds } from './scan.mjs';
export class MoleculeViewer {
  constructor(element,onSelect) {
    this.element=element; this.onSelect=onSelect; this.model=null; this.labels=true;
    this.pair=[]; this.moving=[];
    try {
      if (!globalThis.$3Dmol) throw new Error('Viewer library could not be loaded.');
      this.viewer=globalThis.$3Dmol.createViewer(element,{backgroundColor:'#f0ece1',antialias:true,minimumZoomToDistance:2.6});
      if (!this.viewer) throw new Error('WebGL is unavailable.');
      this.viewer.setProjection('orthographic');
      new ResizeObserver(()=>{
        this.viewer.resize();
        if(this.model)this.fit();
      }).observe(element);
    } catch(error) {
      element.textContent=`3D viewer unavailable: ${error.message} Calculations and atom-index scan controls still work.`;
      element.classList.add('viewer-fallback');
    }
  }
  setMolecule(molecule,{reset=false}={}) {
    this.molecule=molecule;
    if (!this.viewer) return;
    const oldView=this.model&&!reset?this.viewer.getView():null;
    this.viewer.removeAllModels(); this.viewer.removeAllLabels();
    const bonds=inferBonds(molecule);
    const atoms=molecule.atoms.map((atom,index)=>({elem:atom.symbol,index,serial:index+1,
      x:atom.xyz[0]*ANGSTROM_PER_BOHR,y:atom.xyz[1]*ANGSTROM_PER_BOHR,z:atom.xyz[2]*ANGSTROM_PER_BOHR,
      bonds:[],bondOrder:[]}));
    for (const [a,b] of bonds) { atoms[a].bonds.push(b);atoms[a].bondOrder.push(1);atoms[b].bonds.push(a);atoms[b].bondOrder.push(1); }
    this.model=this.viewer.addModel(); this.model.addAtoms(atoms);
    this.viewer.setClickable({},true,atom=>this.onSelect(atom.index));
    this.style();
    if (oldView) this.viewer.setView(oldView); else {
      if (!this.oriented) { this.viewer.rotate(65,'y');this.viewer.rotate(12,'x');this.oriented=true; }
      this.fit();
    }
    this.viewer.render();
    this.element.dataset.atomCount=String(atoms.length);
  }
  style() {
    if (!this.viewer||!this.model) return;
    const style={stick:{radius:.12,colorscheme:'Jmol'},sphere:{scale:.32,colorscheme:'Jmol'}};
    this.viewer.setStyle({},style);
    if (this.moving.length) this.viewer.setStyle({index:this.moving},{stick:{radius:.12,colorscheme:'Jmol'},sphere:{scale:.32,colorscheme:'Jmol'}});
    if (this.pair.length) this.viewer.setStyle({index:this.pair},{stick:{radius:.13,color:'#c45c26'},sphere:{scale:.35,color:'#c45c26'}});
    this.viewer.removeAllLabels();
    if (this.labels) this.molecule.atoms.forEach((atom,i)=>{
      this.viewer.addLabel(`${i+1} ${atom.symbol}`,{position:{x:atom.xyz[0]*ANGSTROM_PER_BOHR,y:atom.xyz[1]*ANGSTROM_PER_BOHR,z:atom.xyz[2]*ANGSTROM_PER_BOHR},
        fontSize:11,fontColor:'#294e59',backgroundColor:'#ffffff',backgroundOpacity:.92,borderRadius:5,inFront:true,screenOffset:{x:13,y:13}});
    });
    this.viewer.render();
  }
  select(pair,moving=[]) { this.pair=pair;this.moving=moving;this.style(); }
  setLabels(value) { this.labels=value;this.style(); }
  fit() {
    if(!this.viewer||!this.model)return;
    this.viewer.zoomTo();
    // Frame the projected structure, with room for spheres, at any canvas aspect ratio.
    const corners=[];
    for(const atom of this.model.selectedAtoms({}))for(const dx of [-.65,.65])for(const dy of [-.65,.65])for(const dz of [-.65,.65])
      corners.push({x:atom.x+dx,y:atom.y+dy,z:atom.z+dz});
    const points=this.viewer.modelToScreen(corners);
    if(points.length){
      const width=Math.max(...points.map(p=>p.x))-Math.min(...points.map(p=>p.x));
      const height=Math.max(...points.map(p=>p.y))-Math.min(...points.map(p=>p.y));
      const scale=Math.min(this.element.clientWidth*.72/width,this.element.clientHeight*.76/height);
      if(Number.isFinite(scale)&&scale>0)this.viewer.zoom(Math.min(8,Math.max(.1,scale)));
    }
    this.viewer.render();
  }
  reset() { this.fit(); }
  png() { return this.viewer?.pngURI(); }
}


