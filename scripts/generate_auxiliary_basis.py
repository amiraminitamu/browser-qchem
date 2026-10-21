"""Regenerate the H-Ne cc-pVDZ-RI correlation-fitting catalog (PySCF 2.14.0)."""
import json
from pathlib import Path
import pyscf
from pyscf import gto
if pyscf.__version__ != '2.14.0':
    raise RuntimeError('Use PySCF 2.14.0 for the pinned auxiliary basis.')
elements = {}
for symbol in ['H','He','Li','Be','B','C','N','O','F','Ne']:
    shells = []
    for angular, *primitives in gto.basis.load('cc-pvdz-ri', symbol):
        for contraction in range(1, len(primitives[0])):
            shells.append({'l':angular, 'exponents':[p[0] for p in primitives],
                           'coefficients':[p[contraction] for p in primitives]})
    elements[symbol] = {'shells':shells}
source = '// SPDX-License-Identifier: AGPL-3.0-only\n'
source += '// cc-pVDZ-RI correlation-fitting parameters from PySCF 2.14.0.\n'
source += '// Spherical functions; general contractions split in PySCF AO order.\n'
source += 'export const AUXILIARY_BASIS='+json.dumps({'id':'cc-pvdz-ri','label':'cc-pVDZ-RI','elements':elements},separators=(',',':'))+';\n'
source += '''export function packAuxiliary(molecule) {
  const shells=[],exponents=[],coefficients=[];let nao=0;
  molecule.atoms.forEach((atom,index)=>{
    const element=AUXILIARY_BASIS.elements[atom.symbol];
    if(!element)throw new Error('No auxiliary basis for '+atom.symbol+'.');
    for(const shell of element.shells){shells.push(index,shell.l,shell.exponents.length);exponents.push(...shell.exponents);coefficients.push(...shell.coefficients);nao+=2*shell.l+1;}
  });
  if(nao>1024)throw new Error('The auxiliary basis exceeds the 1024-function limit.');
  return {shells,exponents,coefficients,nao};
}
'''
(Path(__file__).resolve().parents[1]/'web/auxiliary-basis.mjs').write_text(source)
