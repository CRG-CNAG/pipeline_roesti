#!/usr/bin/env python3
import argparse
import requests
import json
import sys

endpoint = "http://dbspipes.crg.es/api/v1/analyses/"
HEADERS = {'Content-type': 'application/json'}


def update_analysis(analysisId, status, verbose=1):

    # sending get request and saving the response as response object
    URL = endpoint + str(analysisId)

    PARAMS = json.dumps({'status': status})

    try:
        print("Sending update analysis status request, URL={} analysisId={} and status={:d}".format(URL, analysisId, status))
        r = requests.put(url = URL, data=PARAMS, headers=HEADERS)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:  # This is the correct syntax
        print(e)
        sys.exit(1)

        
def update_progress(analysisId, job, n):

    # sending get request and saving the response as response object
    URL = endpoint + str(analysisId) + '/progress'

    PARAMS = json.dumps({"progress": {'job_name': job, 'task_n': n}})

    try:
        print("Sending new progress request, URL={} and progress={} with number={:d}".format(URL, job, n))
        r = requests.post(url = URL, data=PARAMS, headers=HEADERS)
        r.raise_for_status()
    except requests.exceptions.RequestException as e:  # This is the correct syntax
        print(e)
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysisId', default=None, type=str)
    parser.add_argument('--sendMessageToWebServer', action='store_true')
    parser.add_argument('--status', default=-99, type=int)
    parser.add_argument('--type', default=0, type=int)
    parser.add_argument('--progress', default=None)
    parser.add_argument('--n', default=-99, type=int)
    options = parser.parse_args()
    if options.sendMessageToWebServer:
        if options.type == 1 and options.progress is not None:
            update_progress(options.analysisId, options.progress, options.n)
        else:
            update_analysis(options.analysisId, options.status)
