#!/usr/bin/env python3
from socketIO_client import SocketIO, LoggingNamespace
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('--analysisId', default=None, type=str)
parser.add_argument('--sendMessageToWebServer', action='store_true')
options = parser.parse_args()

if options.sendMessageToWebServer:
    with SocketIO('dbspipes.crg.es', 50001, LoggingNamespace) as socketIO:
        # Send message using web socket to the web server DBSpipes
        status = 1
        print("Sending message to web server via websocket, analysisId={:s} and status={:d}".format(options.analysisId, status))
        data = {"internal_id": options.analysisId, "status": status}
        socketIO.emit('on_update', json.dumps(data))
        # Listen
        socketIO.wait(seconds=1)
