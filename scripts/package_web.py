"""Package one cache-consistent browser build after the WASM target is linked."""
import hashlib,re,shutil,sys
from pathlib import Path

def package(web,dist):
    web,dist=Path(web),Path(dist)
    assets=sorted(p for p in web.iterdir() if p.suffix in {'.html','.mjs','.css'})
    engine=[dist/'qc_core.mjs',dist/'qc_core.wasm']
    vendor=sorted(p for p in (web/'vendor').rglob('*') if p.is_file()) if (web/'vendor').exists() else []
    digest=hashlib.sha256()
    for p in [*assets,*vendor,*engine]:
        raw=p.read_bytes()
        name=str(p.relative_to(web)) if web in p.parents else 'engine/'+p.name
        digest.update(name.encode()+b'\0'+str(len(raw)).encode()+b'\0'+raw)
    build=digest.hexdigest()[:20]
    # Rewrite local module imports, Worker URLs, scripts and styles together.
    # Generated qc_core.mjs is versioned by its import; the worker's locateFile
    # callback applies this same build ID to the generated WASM fetch.
    local=re.compile(r"""(['"])(\./[^'"?\s]+\.(?:mjs|js|css|wasm))(?:\?[^'"]*)?\1""")
    for p in assets:
        text=p.read_text()
        text=local.sub(lambda m:m[1]+m[2]+'?v='+build+m[1],text)
        if p.name=='build-info.mjs':
            text=re.sub(r"export const BUILD_ID='[^']*';","export const BUILD_ID='"+build+"';",text)
        (dist/p.name).write_text(text)
    if (web/'vendor').exists():
        shutil.copytree(web/'vendor',dist/'vendor',dirs_exist_ok=True)
    (dist/'build.json').write_text('{"buildId":"'+build+'","protocol":2}\n')
    print('Browser build '+build)
    return build

if __name__=='__main__':
    package(*sys.argv[1:])
