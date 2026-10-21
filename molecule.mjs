// SPDX-License-Identifier: AGPL-3.0-only
import { STO3G, BASIS_CATALOG, MAX_AO } from './basis.mjs?v=e933e699b7cf09943ce2';
export const ANGSTROM_PER_BOHR = 0.52917721092; // PySCF reference conversion
export const WATER_XYZ = `3
Water, coordinates in angstrom
O  0.000000  0.000000  0.000000
H  0.000000 -0.757000  0.587000
H  0.000000  0.757000  0.587000`;

export function parseGeometry(text, unit = 'angstrom', basis = 'sto-3g') {
  basis = basis.toLowerCase();
  if (!BASIS_CATALOG[basis]) throw new Error('Choose a supported basis set.');
  if (!['angstrom', 'bohr'].includes(unit)) throw new Error('Choose angstrom or bohr.');
  let lines = text.replace(/^\uFEFF/, '').trim().split(/\r?\n/);
  let expected;
  if (/^\d+$/.test(lines[0].trim())) {
    expected = Number(lines[0]);
    lines = lines.slice(2); // XYZ comment line is required, possibly blank.
  }
  lines = lines.filter(line => line.trim());
  if (!lines.length || lines.length > 32 || (expected !== undefined && expected !== lines.length)) {
    throw new Error('Enter 1 to 32 atoms; the XYZ count must match the coordinate lines.');
  }
  const atoms = lines.map((line, i) => {
    const parts = line.trim().split(/\s+/);
    if (parts.length !== 4 || !/^[A-Za-z]{1,2}$/.test(parts[0])) throw new Error(`Atom ${i+1}: expected Element x y z.`);
    const symbol = parts[0][0].toUpperCase() + parts[0].slice(1).toLowerCase();
    if (!STO3G[symbol]) throw new Error(`Unsupported element ${symbol}. Supported: H through Ne.`);
    const xyz = parts.slice(1).map(value => {
      if (!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?$/.test(value)) throw new Error(`Atom ${i+1}: coordinates must be finite numbers.`);
      return Number(value.replace(/[dD]/, 'e')) / (unit === 'angstrom' ? ANGSTROM_PER_BOHR : 1);
    });
    if (xyz.some(value => !Number.isFinite(value) || Math.abs(value) > 1e6)) throw new Error(`Atom ${i+1}: invalid or excessively large coordinate.`);
    return { symbol, charge: STO3G[symbol].charge, xyz };
  });
  return { atoms, basis, coordinateUnit: 'bohr' };
}

export function packBasis(molecule, { enforceLimit = true } = {}) {
  const catalog = BASIS_CATALOG[molecule.basis ?? 'sto-3g'];
  if (!catalog) throw new Error('Unsupported basis set.');
  const charges = [], xyz = [], shells = [], exponents = [], coefficients = [], aoLabels = [];
  molecule.atoms.forEach((atom, index) => {
    const entry = catalog.elements[atom.symbol];
    if (!entry) throw new Error(`Unsupported element ${atom.symbol}.`);
    charges.push(entry.charge); xyz.push(...atom.xyz);
    const shellCounts = [0,0,0];
    entry.shells.forEach(shell => {
      shells.push(index, shell.l, shell.exponents.length);
      exponents.push(...shell.exponents); coefficients.push(...shell.coefficients);
      const count = ++shellCounts[shell.l];
      const labels = shell.l === 0 ? [`s${count}`] : shell.l === 1
        ? ['x','y','z'].map(axis => `p${count}${axis}`)
        : ['xy','yz','z2','xz','x2-y2'].map(component => `d${count}${component}`);
      aoLabels.push(...labels.map(label => `${index+1} ${atom.symbol} ${label}`));
    });
  });
  if (enforceLimit && aoLabels.length > MAX_AO) throw new Error(`This basis needs ${aoLabels.length} AOs; the current limit is ${MAX_AO} AOs. Use a smaller basis or molecule.`);
  return { charges, xyz, shells, exponents, coefficients, aoLabels };
}

export function basisSummary(molecule) {
  const data = packBasis(molecule,{ enforceLimit:false });
  const nao = data.aoLabels.length;
  return { nao, nshell:data.shells.length/3, eriBytes:8*nao**4, fits:nao<=MAX_AO, limit:MAX_AO };
}

export function toXYZ(molecule, unit = 'angstrom', comment = '') {
  const scale = unit === 'angstrom' ? ANGSTROM_PER_BOHR : 1;
  return [molecule.atoms.length,String(comment).replace(/[\r\n]/g,' '),
    ...molecule.atoms.map(atom => `${atom.symbol.padEnd(2)} ${atom.xyz.map(x => (x*scale).toFixed(12)).join('  ')}`)].join('\n')+'\n';
}
