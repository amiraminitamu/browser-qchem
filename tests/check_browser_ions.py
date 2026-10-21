"""Check ionic presets, spin/charge validation, UHF scans and exports in Chromium."""
from browser_controls import control

import argparse,csv,functools,http.server,json,threading
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
fixture=json.loads((ROOT/'tests/fixtures/methods.json').read_text())
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--browser-path');args=parser.parse_args()
class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self,*args):pass
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Quiet,directory=str(ROOT/'build/wasm/dist')))
threading.Thread(target=server.serve_forever,daemon=True).start()
try:
 with sync_playwright() as p:
    launch={'headless':True,'args':['--enable-unsafe-swiftshader']}
    if args.browser_path:launch['executable_path']=args.browser_path
    browser=p.chromium.launch(**launch);page=browser.new_page(viewport={'width':1512,'height':982},reduced_motion='reduce');errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(f'http://127.0.0.1:{server.server_port}')
    def status(prefix):page.wait_for_function('(p)=>document.querySelector("#status").textContent.startsWith(p)',arg=prefix,polling=50)
    def download(selector):
        with page.expect_download() as d:control(page, selector).click()
        return Path(d.value.path()).read_text()
    def reference(name,basis,method='uhf'):return next(c for c in fixture['cases'] if c['name']==name and c['basis']==basis and c['options']['method']==method)
    status('CONVERGED')
    control(page, '#scf-integrals').select_option('exact')
    for preset,name,charge,mult,ne in [('h2plus','H2+',1,2,1),('heh','HeH',0,2,3),('hehplus','HeH+',1,1,2)]:
        control(page, '#preset').select_option(preset)
        assert control(page, '#charge').input_value()==str(charge)
        assert control(page, '#multiplicity').input_value()==str(mult)
        if preset!='hehplus':assert page.locator('#method').input_value()=='uhf'
        else:
            assert page.locator('#method').input_value()=='rhf'
            page.keyboard.press('Control+Enter');status('CONVERGED')
            closed=json.loads(download('#download'));assert closed['rhf']['charge']==1 and closed['rhf']['electrons']==2
            page.locator('#method').select_option('uhf')
        for basis in ['sto-3g','6-31g','cc-pvdz']:
            page.locator('#basis').select_option(basis);page.keyboard.press('Control+Enter');status('CONVERGED')
            result=json.loads(download('#download'));u=result['uhf'];ref=reference(name,basis)
            assert u['electrons']==ne and u['charge']==charge and u['multiplicity']==mult
            assert abs(u['energy']-ref['referenceEnergy'])<2e-8
            assert u['nalpha']+u['nbeta']==ne
            assert abs(u['spinSquare']-ref['spinSquare'])<2e-6
            assert page.locator('#metric-electrons').inner_text()==str(ne)
            if charge:assert '⁺' in page.locator('#formula').inner_text()
    # Load the exact one-electron scan geometry, with parameters in bohr.
    control(page, '#preset').select_option('h2plus');page.locator('#basis').select_option('6-31g');control(page, '#mode').select_option('scan')
    for key,value in [('a','1'),('b','2'),('moving','2'),('start','1.4'),('end','4'),('count','3')]:control(page, '#scan-'+key).fill(value)
    control(page, '#scan-unit').select_option('bohr');page.keyboard.press('Control+Enter')
    page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete")',polling=50)
    scan=json.loads(download('#scan-json'));ref=next(c for c in fixture['ionScans'] if c['name']=='H2+' and c['basis']=='6-31g')
    assert scan['spec']['charge']==1 and scan['spec']['electrons']==1
    for pt,energy in zip(scan['points'],ref['energies']):assert abs(pt['energy']-energy)<2e-8
    rows=list(csv.DictReader(download('#scan-csv').splitlines()));assert all(r['charge']=='1' and r['electrons']=='1' and r['multiplicity']=='2' for r in rows)
    assert 'charge=1 electrons=1' in download('#scan-xyz')
    page.locator('#scan-point-table tbody tr').nth(1).click();assert control(page, '#charge').input_value()=='1' and control(page, '#multiplicity').input_value()=='2'
    control(page, '#mode').select_option('rhf');page.locator('#method').select_option('mp2');page.keyboard.press('Control+Enter');status('Input or calculation error');assert 'even electron count' in page.locator('#summary').inner_text()
    page.locator('#method').select_option('uhf');control(page, '#multiplicity').fill('1');page.keyboard.press('Control+Enter');status('Input or calculation error');assert 'Multiplicity' in page.locator('#summary').inner_text()
    # Test an anion and verify charge superscripts update independently of preset.
    page.locator('#coordinates-editor').evaluate('(e)=>e.open=true');control(page, '#geometry').fill('H 0 0 0');control(page, '#charge').fill('-1');page.locator('#basis').select_option('6-31g');page.keyboard.press('Control+Enter');status('CONVERGED')
    result=json.loads(download('#download'));assert result['uhf']['electrons']==2 and result['uhf']['charge']==-1
    assert abs(result['uhf']['energy']-reference('H-','6-31g')['referenceEnergy'])<2e-8
    assert page.locator('#formula').inner_text()=='H⁻'
    page.set_viewport_size({'width':390,'height':844});page.locator('#molecule-popover > summary').click();control(page, '#preset').select_option('heh');page.locator('#calculate').click();status('CONVERGED')
    assert page.locator('#metric-electrons').inner_text()=='3';assert 'doublet' in page.locator('#molecule-state').inner_text()
    assert not errors,errors;browser.close();print('H2+/HeH/HeH+ presets and UHF energies in three bases, charged PES/CSV/XYZ/JSON, anions, spin rejection and mobile recovery PASS')
finally:server.shutdown();server.server_close()
