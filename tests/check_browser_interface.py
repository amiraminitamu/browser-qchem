"""Real-WASM checks for the redesigned navigation, controls and responsive interface."""
import argparse,functools,http.server,threading,io
from PIL import Image
from pathlib import Path
from playwright.sync_api import sync_playwright
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--browser-path');args=parser.parse_args()
class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self,*args):pass
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),functools.partial(Quiet,directory=str(ROOT/'build/wasm/dist')))
threading.Thread(target=server.serve_forever,daemon=True).start()
try:
    with sync_playwright() as p:
        launch={'headless':True,'args':['--enable-unsafe-swiftshader']}
        if args.browser_path:launch['executable_path']=args.browser_path
        browser=p.chromium.launch(**launch)
        page=browser.new_page(viewport={'width':1512,'height':900},reduced_motion='reduce')
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(f'http://127.0.0.1:{server.server_port}')
        def status(prefix):page.wait_for_function('(p)=>document.querySelector("#status").textContent.startsWith(p)',arg=prefix,polling=50)
        status('CONVERGED')
        assert page.locator('#job-list .job-item').count()>=1
        assert page.locator('.menu-bar').is_visible()
        assert page.locator('.app-toolbar').is_visible()
        assert page.locator('.brand-name').inner_text()=='Lorentzian'
        page.locator('.menu-bar .menu-button', has_text='Analysis').click()
        page.locator('[data-menu="tab"][data-tab-name="orbitals"]').click()
        assert page.locator('#tab-orbitals').get_attribute('aria-selected')=='true'
        page.locator('#tab-scf').click()
        assert not page.locator('#coordinates-editor').evaluate('(e)=>e.open')
        page.locator('#nav-molecule').click();assert page.locator('#coordinates-editor').evaluate('(e)=>e.open')
        page.locator('#coordinates-editor summary').click()
        assert not page.locator('.advanced-settings').evaluate('(e)=>e.open')
        assert not page.locator('.summary-details').evaluate('(e)=>e.open')
        run_y=page.locator('#calculate').bounding_box()['y']
        assert 0<run_y<850,'Primary action should stay within the viewport'
        assert page.locator('#metric-aos').inner_text()=='7'
        viewer_box=page.locator('.viewer-card').bounding_box();energy_box=page.locator('.energy-card').bounding_box()
        assert energy_box['x']>=viewer_box['x']+viewer_box['width'],'Structure and energy should sit side by side'
        assert page.locator('#metric-electrons').inner_text()=='10'
        page.locator('#mode').select_option('scan')
        assert page.locator('#mode').input_value()=='scan'
        assert page.locator('#scan-a').is_visible()
        assert abs(page.locator('#calculate').bounding_box()['y']-run_y)<1
        assert page.locator('#run-label').inner_text()=='Run bond scan'
        page.locator('#mode').select_option('rhf')
        page.locator('.advanced-settings summary').click()
        page.locator('#maxiter').fill('1');page.keyboard.press('Control+Enter');status('NOT CONVERGED')
        assert page.locator('.summary-details').evaluate('(e)=>e.open')
        assert 'Unconverged energy' in page.locator('#summary').inner_text()
        page.locator('#maxiter').fill('100');page.keyboard.press('Control+Enter');status('CONVERGED')
        page.locator('#scf-chart circle[role=button]').first.click()
        assert page.locator('#iteration-details').evaluate('(e)=>e.open')
        assert page.locator('#scf-history tr.selected').count()==1
        page.locator('#tab-scf').click();page.keyboard.press('ArrowRight')
        assert page.locator('#tab-orbitals').get_attribute('aria-selected')=='true'
        page.keyboard.press('End');assert page.locator('#tab-data').get_attribute('aria-selected')=='true'
        page.keyboard.press('Home');assert page.locator('#tab-scf').get_attribute('aria-selected')=='true'
        page.locator('#help-open').click();assert page.locator('#workspace-guide').evaluate('(e)=>e.open')
        page.keyboard.press('Escape');assert not page.locator('#workspace-guide').evaluate('(e)=>e.open')
        page.evaluate('scrollTo(0,0)');page.locator('#focus-view').click()
        assert page.locator('#focus-view').get_attribute('aria-pressed')=='true'
        assert page.locator('#molecule-viewer').bounding_box()['height']>500
        page.keyboard.press('Escape');assert page.locator('#focus-view').get_attribute('aria-pressed')=='false'
        page.evaluate('scrollTo(0,0)')
        page.evaluate("document.querySelector('.setup-scroll').scrollTop=0;document.querySelector('.advanced-settings').open=false;document.querySelector('#coordinates-editor').open=false;document.querySelector('#iteration-details').open=false;document.activeElement.blur()")
        page.screenshot(path=str(ROOT.parent/'redesign-desktop.png'),full_page=True)
        print('Scientific layout, persistent Run, scan setup, convergence settings, keyboard run, tabs, diagnostics, help and expanded viewer PASS',flush=True)
        page.set_viewport_size({'width':390,'height':844});page.evaluate('scrollTo(0,0)')
        assert page.locator('#setup-panel').get_attribute('aria-hidden')=='true'
        assert page.locator('#mobile-run').is_visible()
        page.locator('#mobile-setup-toggle').click()
        assert page.locator('#setup-panel').bounding_box()['x']>=0
        assert page.locator('#mobile-setup-toggle').get_attribute('aria-expanded')=='true'
        page.screenshot(path=str(ROOT.parent/'redesign-drawer.png'))
        page.locator('#preset').select_option('h2');page.locator('#basis').select_option('6-31g')
        page.locator('#calculate').click();status('CONVERGED')
        assert page.locator('#setup-panel').get_attribute('aria-hidden')=='true'
        assert page.locator('#metric-aos').inner_text()=='4'
        page.locator('#mobile-edit').click();page.locator('#mode').select_option('scan')
        for key,value in [('a','1'),('b','2'),('moving','2'),('start','.5'),('end','1.5'),('count','3')]:page.locator('#scan-'+key).fill(value)
        page.locator('#calculate').click()
        page.wait_for_function('document.querySelector("#progress-text").textContent.startsWith("Complete")',polling=50)
        assert page.locator('#scan-point-table tbody tr').count()==3
        page.locator('#pes-chart circle[role=button]').nth(1).click()
        assert '1.000000000000' in page.locator('#geometry').input_value()
        assert abs(float(page.locator('#energy-value').inner_text())-(-1.094807962861))<1e-9
        image=Image.open(io.BytesIO(page.locator('#molecule-viewer canvas').screenshot())).convert('RGB')
        visible=[(x,y) for y in range(image.height) for x in range(image.width) if (lambda rgb:max(rgb)>90 and max(rgb)-min(rgb)>45)(image.getpixel((x,y)))]
        assert visible
        assert max(x for x,y in visible)-min(x for x,y in visible)>70,'Molecule should remain clearly visible on a phone'
        page.evaluate('scrollTo(0,0)')
        page.screenshot(path=str(ROOT.parent/'redesign-mobile.png'),full_page=True)
        for width in [360,390,768,1024,1440]:
            page.set_viewport_size({'width':width,'height':900})
            assert page.evaluate('document.documentElement.scrollWidth<=innerWidth'),width
            assert page.locator('#energy-value').bounding_box()['height']<60,width
        assert not errors,errors
        browser.close()
        print('Mobile drawer, real single-point/scan calculations, curve inspection and responsive widths 360-1440px PASS')
finally:
    server.shutdown();server.server_close()

