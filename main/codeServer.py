# Python 3 server example
import os
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import parse_qs
from subprocess import PIPE, STDOUT, run

hostName = "localhost"
serverPort = 8080

class QuizRequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        self.wfile.write(bytes(read_template("index.html"), "utf-8"))

    def do_POST(self):
        data_string = self.rfile.read(int(self.headers['Content-Length'])).decode('utf-8')
        data = parse_qs(data_string, keep_blank_values=True)

        self.send_response(200)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        if self.path == "/run_code":
            code = data['codestuff'][0]
            p = run("python3", stdout=PIPE, shell=True, stderr=STDOUT, input=code, encoding='ascii')
            output = p.stdout

            new_page = read_template("index.html").replace("<!-- OUTPUT PLACEHOLDER -->", output, 1)
            new_page = new_page.replace("ENTER CODE HERE", code)
            self.wfile.write(bytes(new_page, "utf-8"))


def read_template(filename, directory='templates'):
    pathname = os.path.join(directory, filename)
    f = open(pathname, "r", encoding="utf-8")
    return f.read()


if __name__ == "__main__":
    webServer = HTTPServer((hostName, serverPort), QuizRequestHandler)
    print("Server has now started up http://%s:%s" % (hostName, serverPort))

    try:
        webServer.serve_forever()
    except KeyboardInterrupt:
        pass

    webServer.server_close()
    print("Server terminated.")
