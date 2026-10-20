"""Serve the compiled workbench with revalidation disabled for all assets.
python3 scripts/serve.py --port 8000 --directory build/wasm/dist
"""
import argparse,functools,http.server

class NoCacheHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cache-Control','no-store, no-cache, must-revalidate')
        self.send_header('Pragma','no-cache')
        self.send_header('Expires','0')
        super().end_headers()

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port',type=int,default=8000)
    parser.add_argument('--directory',default='build/wasm/dist')
    args=parser.parse_args()
    server=http.server.ThreadingHTTPServer(('127.0.0.1',args.port),functools.partial(NoCacheHandler,directory=args.directory))
    print('Open http://localhost:'+str(args.port),flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:server.server_close()
