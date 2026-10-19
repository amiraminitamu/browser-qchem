// SPDX-License-Identifier: AGPL-3.0-only
import { STO3G } from './basis.mjs';
export const ANGSTROM_PER_BOHR = 0.52917721092; // PySCF reference conversion
export const WATER_XYZ = `3
Water, coordinates in angstrom
O  0.000000  0.000000  0.000000
H  0.000000 -0.757000  0.587000
H  0.000000  0.757000  0.587000`;

export function parseGeometry(text, unit = 'angstrom') {
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
    if (!STO3G[symbol]) throw new Error(`STO-3G data for ${symbol} are not included. Supported: H through Ne.`);
    const xyz = parts.slice(1).map(value => {
      if (!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?$/.test(value)) throw new Error(`Atom ${i+1}: coordinates must be finite numbers.`);
      return Number(value.replace(/[dD]/, 'e')) / (unit === 'angstrom' ? ANGSTROM_PER_BOHR : 1);
    });
    if (xyz.some(value => !Number.isFinite(value) || Math.abs(value) > 1e6)) throw new Error(`Atom ${i+1}: invalid or excessively large coordinate.`);
    return { symbol, charge: STO3G[symbol].charge, xyz };
  });
  return { atoms, basis: 'sto-3g', coordinateUnit: 'bohr' };
}

export function packBasis(molecule) {
  const charges = [], xyz = [], shells = [], exponents = [], coefficients = [], aoLabels = [];
  molecule.atoms.forEach((atom, index) => {
    const entry = STO3G[atom.symbol];
    if (!entry) throw new Error(`Unsupported element ${atom.symbol}.`);
    charges.push(entry.charge); xyz.push(...atom.xyz);
    let sShell = 0;
    entry.shells.forEach(shell => {
      shells.push(index, shell.l, shell.exponents.length);
      exponents.push(...shell.exponents); coefficients.push(...shell.coefficients);
      const labels = shell.l === 0 ? [`${++sShell}s`] : ['2px', '2py', '2pz'];
      aoLabels.push(...labels.map(label => `${index+1} ${atom.symbol} ${label}`));
    });
  });
  if (aoLabels.length > 40) throw new Error('This dense ERI prototype is limited to 40 AOs.');
  return { charges, xyz, shells, exponents, coefficients, aoLabels };
}
