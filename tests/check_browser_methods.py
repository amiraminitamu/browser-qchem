"""Real-WASM browser flows for UHF and conventional RHF-MP2."""
from browser_controls import control

import argparse,functools,http.server,threading,json
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
MP2_WATER=next(c['energy'] for c in json.loads((ROOT/'tests/fixtures/methods.json').read_text())['cases'] if c['name']=='water' and c['basis']=='sto-3g' and c['options']['method']=='mp2')
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
    def run():page.keyboard.press('Control+Enter')
    def geometry(text):
        page.locator('#coordinates-editor').evaluate('(e)=>e.open=true');control(page, '#geometry').fill(text)
    def export():
        with page.expect_download() as d:control(page, '#download').click()
        return json.loads(Path(d.value.path()).read_text())
    status('CONVERGED');page.locator('#method').select_option('uhf');geometry('H 0 0 0');control(page, '#multiplicity').fill('2');run();status('CONVERGED')
    assert 'UHF' in page.locator('#energy-label').inner_text();assert '0.750000' in control(page, '#method-details').inner_text()
    result=export();assert result['uhf']['nalpha']==1 and result['uhf']['nbeta']==0;assert isinstance(result['uhf']['densityAlpha'],list)
    page.locator('#tab-orbitals').click();assert 'Alpha' in page.locator('#orbital-table').inner_text();assert 'Beta' in page.locator('#orbital-table').inner_text()
    page.locator('#tab-matrices').click();page.locator('#matrix-kind').select_option('fockBeta');assert page.locator('#matrix table').count()==1
    control(page, '#multiplicity').fill('1');run();status('Input or calculation error');assert 'Multiplicity' in page.locator('#summary').inner_text()
    geometry('H 0 0 0\nH 0 0 3');run();status('CONVERGED');assert float(page.locator('#energy-value').inner_text())<-.93
    assert export()['uhf']['spinSquare']>.9
    page.locator('#method').select_option('mp2');control(page, '#preset').select_option('water');run();status('CONVERGED')
    result=export();assert abs(result['mp2']['energy']-MP2_WATER)<2e-9
    assert 'MP2' in page.locator('#energy-label').inner_text();assert 'RHF reference' in control(page, '#method-details').inner_text();assert 'MP2 correlation' in control(page, '#method-details').inner_text()
    page.locator('#tab-scf').click();page.evaluate("scrollTo(0,0);document.querySelector('#coordinates-editor').open=false;document.activeElement.blur()");page.screenshot(path=str(ROOT.parent/'methods-workspace.png'))
    page.locator('#tab-orbitals').click();assert 'RHF reference' in page.locator('#orbital-note').inner_text()
    page.locator('.advanced-settings').evaluate('(e)=>e.open=true');control(page, '#maxiter').fill('1');run();status('NOT CONVERGED');result=export();assert result['mp2']['energy'] is None;assert 'not computed' in page.locator('#status').inner_text()
    control(page, '#maxiter').fill('100');control(page, '#preset').select_option('h2');control(page, '#mode').select_option('scan')
    for key,value in [('a','1'),('b','2'),('moving','2'),('start','.7'),('end','1.3'),('count','3')]:control(page, '#scan-'+key).fill(value)
    run();page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete")',polling=50)
    assert page.locator('#scan-point-table tbody tr').count()==3;assert 'MP2' in page.locator('#pes-note').inner_text()
    with page.expect_download() as d:control(page, '#scan-json').click()
    scan=json.loads(Path(d.value.path()).read_text());assert all(abs(pt['energy']-pt['result']['mp2']['energy'])<1e-12 for pt in scan['points'])
    page.locator('#scan-point-table tbody tr').nth(1).click();assert page.locator('#method').input_value()=='mp2';assert 'MP2' in page.locator('#energy-label').inner_text()
    page.locator('#method').select_option('uhf');control(page, '#multiplicity').fill('3');run()
    page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete")',polling=50)
    with page.expect_download() as d:control(page, '#scan-json').click()
    scan=json.loads(Path(d.value.path()).read_text());assert scan['spec']['method']=='uhf';assert scan['spec']['multiplicity']==3
    assert all(pt['result']['uhf']['nalpha']==2 and pt['result']['uhf']['nbeta']==0 for pt in scan['points'])
    page.locator('#method').select_option('rhf');page.locator('#scan-point-table tbody tr').first.click()
    assert page.locator('#method').input_value()=='rhf' and control(page, '#multiplicity').input_value()=='3'
    assert 'current setup RHF/' in page.locator('#status').inner_text()
    page.screenshot(path=str(ROOT.parent/'methods-desktop.png'),full_page=True)
    page.set_viewport_size({'width':390,'height':844});page.locator('#molecule-popover > summary').click();page.locator('#method').select_option('uhf');control(page, '#multiplicity').fill('1');control(page, '#guess').select_option('mix');control(page, '#mode').select_option('rhf');control(page, '#preset').select_option('h2');page.locator('#calculate').click();status('CONVERGED');assert 'UHF' in page.locator('#energy-label').inner_text()
    page.screenshot(path=str(ROOT.parent/'methods-mobile.png'),full_page=True)
    assert control(page, '#method-details').is_visible();assert page.evaluate('document.documentElement.scrollWidth<=innerWidth');assert not errors,errors
    browser.close();print('UHF spins/arrays/matrices, broken symmetry, MP2 components/reference failure, MP2 scan/export, mobile method controls PASS')
finally:server.shutdown()
