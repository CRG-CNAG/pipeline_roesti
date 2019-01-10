#!/usr/bin/env python3
from socketIO_client import SocketIO, LoggingNamespace
import argparse
import json


def send_socket_message(analysisId, status, verbose=1):
    with SocketIO('dbspipes.crg.es', 50001, LoggingNamespace) as socketIO:
        # Send message using web socket to the web server DBSpipes
        if verbose >= 1:
            print("Sending message to web server via websocket, analysisId={} and status={:d}".format(analysisId, status))
        data = {"internal_id": analysisId, "status": status}
        socketIO.emit('on_update', json.dumps(data))
        # Listen
        socketIO.wait(seconds=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysisId', default=None, type=str)
    parser.add_argument('--sendMessageToWebServer', action='store_true')
    parser.add_argument('--status', default=-99, type=int)
    options = parser.parse_args()
    if options.sendMessageToWebServer:
        send_socket_message(options.analysisId, options.status)
