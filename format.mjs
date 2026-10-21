// SPDX-License-Identifier: AGPL-3.0-only
// Apply to display text only: machine-readable files and input values stay ASCII.
export function mathText(value) {return String(value).replace(/(^|[\s=(eE])-(?=\d)/g,'$1−').replace(/\bEh\b/g,'Eₕ');}
