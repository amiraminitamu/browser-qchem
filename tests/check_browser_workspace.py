"""End-to-end basis, viewer and PES tests in Chromium, using the real WASM worker."""
import argparse,csv,functools,http.server,io,json,math,tempfile,threading
from pathlib import Path
from playwright.sync_api import sync_playwright
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--browser-path');args=parser.parse_args()
class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self,*args): pass
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(QuietHandler,directory=str(ROOT/'build/wasm/dist')))
threading.Thread(target=server.serve_forever,daemon=True).start()
reference=json.loads((ROOT/'tests/fixtures/workspace.json').read_text())['cases']
try:
    with sync_playwright() as p,tempfile.TemporaryDirectory() as temp:
        launch={'headless':True,'args':['--enable-unsafe-swiftshader']}
        if args.browser_path:launch['executable_path']=args.browser_path
        browser=p.chromium.launch(**launch)
        page=browser.new_page(viewport={'width':1512,'height':1100})
        errors=[];external=[];page.on('pageerror',lambda e:(errors.append(str(e)),print('BROWSER ERROR',e,flush=True)))
        page.on('request',lambda r:external.append(r.url) if not r.url.startswith(('http://127.0.0.1:','blob:','data:')) else None)
        # Expose only the renderer wrapper in this test response to project atom centers.
        # The atom selection itself is exercised with actual canvas mouse clicks.
        page.route('**/app.mjs*',lambda route:route.fulfill(response=route.fetch(),body=(ROOT/'build/wasm/dist/app.mjs').read_text().replace('const viewer=new MoleculeViewer','const viewer=window.testViewer=new MoleculeViewer')))
        page.goto(f'http://127.0.0.1:{server.server_port}')
        def status(prefix):
            page.wait_for_function('(p)=>document.querySelector("#status").textContent.startsWith(p)',arg=prefix,timeout=30000,polling=100)
        def idle():
            page.wait_for_function('!document.querySelector("#calculate").disabled',timeout=30000,polling=100)
        def download(selector,filename):
            with page.expect_download() as d:page.locator(selector).click()
            target=Path(temp)/filename;d.value.save_as(target);return target
        def scan_complete():
            page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete") && !document.querySelector("#calculate").disabled',timeout=30000,polling=100)
        status('CONVERGED')
        assert page.locator('#molecule-viewer canvas').count()==1
        assert page.locator('#molecule-viewer').get_attribute('data-atom-count')=='3'
        image=download('#viewer-png','water.png')
        im=Image.open(image).convert('RGB');assert im.width>400 and im.height>200
        colored=sum(1 for y in range(im.height) for x in range(im.width) if (lambda rgb:max(rgb)>90 and max(rgb)-min(rgb)>45)(im.getpixel((x,y))))
        assert colored>1000,'Molecule should occupy a visible area, rather than an empty canvas'
        canvas=page.locator('#molecule-viewer canvas');before=canvas.screenshot()
        box=canvas.bounding_box();page.mouse.move(box['x']+120,box['y']+100);page.mouse.down();page.mouse.move(box['x']+230,box['y']+150,steps=8);page.mouse.up()
        assert canvas.screenshot()!=before,'Dragging should rotate the molecule'
        page.locator('#fit-view').click()
        for index in [2,0]:
            position=page.evaluate('(i)=>testViewer.viewer.modelToScreen(testViewer.model.selectedAtoms({index:i})[0])',index)
            page.mouse.click(position['x'],position['y'])
        assert page.locator('#scan-a').input_value()=='3'
        assert page.locator('#scan-b').input_value()=='1'
        page.locator('#mode').select_option('scan')
        page.locator('#auto-fragment').click()
        assert page.locator('#scan-moving').input_value()=='1, 2'
        print('WebGL rendering, PNG, rotation, atom picking and connected fragment selection PASS',flush=True)
        # The larger basis must be applied to the worker, not just the displayed label.
        page.locator('#mode').select_option('rhf');page.locator('#basis').select_option('cc-pvdz');page.locator('#calculate').click();status('CONVERGED')
        data=json.loads(download('#download','water.json').read_text())
        expected=next(c for c in reference if c['name']=='water' and c['basis']=='cc-pvdz')
        assert data['nao']==24 and abs(data['rhf']['energy']-expected['energy'])<1e-9
        page.locator('#nav-molecule').click()
        page.locator('#geometry').fill('Ne 0 0 0\nNe 0 0 3\nNe 0 0 6')
        page.wait_for_function('document.querySelector("#input-message").textContent.includes("40 AO")',polling=100)
        assert page.locator('#calculate').is_disabled()
        # Upload a real XYZ file, then calculate a five-point H2 curve.
        page.locator('#xyz-file').set_input_files({'name':'h2.xyz','mimeType':'text/plain','buffer':b'2\nHydrogen\nH 0 0 0\nH 0 0 .74\n'})
        page.wait_for_function('document.querySelector("#molecule-viewer").dataset.atomCount==="2"',polling=100)
        page.locator('#basis').select_option('6-31g');page.locator('#mode').select_option('scan')
        for key,value in [('a','1'),('b','2'),('moving','2'),('start','.5'),('end','1.5'),('count','5')]:page.locator('#scan-'+key).fill(value)
        page.locator('#scan-unit').select_option('angstrom');page.locator('#calculate').click();scan_complete()
        assert page.locator('#scan-point-table tbody tr').count()==5
        curve=json.loads(download('#scan-json','scan.json').read_text())
        csvrows=list(csv.DictReader(download('#scan-csv','scan.csv').read_text().splitlines()))
        xyz=download('#scan-xyz','scan.xyz').read_text().strip().splitlines()
        assert len(csvrows)==5 and len(xyz)==20
        for point,row in zip(curve['points'],csvrows):
            expected=next(c for c in reference if c['basis']=='6-31g' and c['name'].startswith('h2-') and float(c['name'][3:])==point['coordinate'])
            assert point['converged'] and abs(point['energy']-expected['energy'])<1e-9
            assert float(row['energy_hartree'])==point['energy'] and row['basis']=='6-31g'
            assert 'eri' not in point['result']
        page.locator('#pes-chart circle[role=button]').nth(1).click()
        assert 'Scan point 2' in page.locator('#results-title').inner_text()
        assert '.750000000000' in page.locator('#geometry').input_value()
        assert abs(float(page.locator('#energy-value').inner_text())-curve['points'][1]['energy'])<1e-9
        assert page.locator('#molecule-viewer').get_attribute('data-atom-count')=='2'
        page.locator('#pes-energy-unit').select_option('absolute-eh')
        page.locator('#scan-point-table tbody tr').nth(3).press('Enter')
        assert 'Scan point 4' in page.locator('#results-title').inner_text()
        assert '1.250000000000' in page.locator('#geometry').input_value()
        page.screenshot(path=str(ROOT.parent/'workspace-scan.png'),full_page=True)
        print('cc-pVDZ worker results, AO preflight, XYZ upload, PES curve, point inspection and all exports PASS',flush=True)
        # Use a larger molecule/basis so completed points remain inspectable before stopping.
        page.locator('#preset').select_option('water');page.locator('#basis').select_option('cc-pvdz')
        page.locator('#scan-start').fill('.8');page.locator('#scan-end').fill('1.3')
        page.locator('#scan-count').fill('101')
        page.locator('#calculate').scroll_into_view_if_needed()
        run_box=page.locator('#calculate').bounding_box()
        run_x=run_box['x']+run_box['width']/2;run_y=run_box['y']+run_box['height']/2
        page.mouse.click(run_x,run_y)
        page.wait_for_function('document.querySelectorAll("#scan-point-table tbody tr").length>=2',polling=50)
        assert page.evaluate('([x,y])=>document.elementFromPoint(x,y).id',[run_x,run_y])=='cancel'
        page.mouse.click(run_x,run_y);idle()
        stopped=json.loads(download('#scan-json','stopped.json').read_text())
        assert stopped['cancelled'] and 1<=len(stopped['points'])<101
        assert 'Stopped' in page.locator('#progress-text').inner_text()
        # Nonconvergence is explicit in both the table and the exported point metadata.
        page.locator('#preset').select_option('water');page.locator('#basis').select_option('sto-3g');page.locator('#scan-a').fill('1');page.locator('#scan-b').fill('2');page.locator('#scan-moving').fill('2')
        page.locator('.advanced-settings').evaluate('(e)=>e.open=true');page.locator('#scan-count').fill('3');page.locator('#scan-start').fill('.8');page.locator('#scan-end').fill('1.2');page.locator('#maxiter').fill('1')
        page.locator('#calculate').click();scan_complete()
        limited=json.loads(download('#scan-json','limited.json').read_text())
        assert all(p['status']=='not_converged' and not p['converged'] for p in limited['points'])
        assert page.locator('#result-badge').inner_text()=='NOT CONVERGED'
        # A colliding intermediate geometry must not leave stale matrices under a failed point.
        page.locator('#geometry').fill('H 0 0 0\nH 1 0 0\nH 2 0 0');page.locator('#charge').fill('1');page.locator('#maxiter').fill('100');page.locator('#scan-start').fill('1');page.locator('#scan-end').fill('3')
        page.locator('#calculate').click();scan_complete()
        page.locator('#scan-point-table tbody tr').nth(1).click()
        assert page.locator('#result-badge').inner_text()=='POINT FAILED'
        assert 'coincide' in page.locator('#summary').inner_text()
        assert page.locator('#matrix tbody tr').count()==0 and page.locator('#download').is_disabled()
        # Disabled hidden scan fields must not block a later single-point form submission.
        page.locator('#scan-start').fill('-1');page.locator('#mode').select_option('rhf');page.locator('#water').click();status('CONVERGED')
        for tab in ['orbitals','matrices','pes','data','scf']:page.locator(f'[data-tab={tab}]').click()
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
        page.screenshot(path=str(ROOT.parent/'workspace-mobile.png'),full_page=True)
        assert not errors,errors
        assert not external,external
        browser.close()
        print('Cancellation, retained exports, nonconvergence, failed-point cleanup, recovery and mobile layout PASS')
        print('No JavaScript exceptions or external asset requests PASS')
finally:
    server.shutdown();server.server_close()

