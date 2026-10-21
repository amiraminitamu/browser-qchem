// SPDX-License-Identifier: AGPL-3.0-only
export const METHOD_SPECS=Object.freeze({
  rhf:{label:'RHF',reference:'RHF',groups:['scf']},
  uhf:{label:'UHF',reference:'UHF',groups:['scf','uhf']},
  mp2:{label:'MP2',reference:'RHF',postSCF:true,groups:['scf','mp2']},
  'df-rmp2':{label:'DF-RMP2',reference:'RHF',postSCF:true,densityFitted:true,groups:['scf','mp2','df']},
  'df-ump2':{label:'DF-UMP2',reference:'UHF',postSCF:true,densityFitted:true,groups:['scf','uhf','mp2','df']}
});
export const isUnrestricted=method=>METHOD_SPECS[method]?.reference==='UHF';
export const isMP2=method=>!!METHOD_SPECS[method]?.postSCF;
