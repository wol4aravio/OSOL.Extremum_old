from flask import Flask, request
import logging

from OSOL_Extremum.computational_core.computational_core import *
from OSOL_Extremum.arithmetics.interval import Interval


def create_app():

    app = Flask(__name__)
    log = logging.getLogger('werkzeug')
    log.disabled = True
    app.logger.disabled = True

    @app.route('/process_request', methods=['GET'])
    def process_request():
        args = json.loads(json.dumps(request.args))
        field = args.pop('field')
        if field == 'target' or field == 'f' or field == 'sim' or field == 'df_grad':
            scope = args.pop('scope', 'real')
            if scope == 'interval':
                for k in args.keys():
                    args[k] = Interval.from_json(args[k])
            elif scope == 'real':
                for k in args.keys():
                    args[k] = float(args[k])
            else:
                raise Exception('Unsupported scope: {}'.format(scope))
            return json.dumps(app.core.request(field, args))
        elif field == 'sim_out':
            json_file = args['json_file']
            save_loc = args['save_loc']
            return json.dumps(app.core.request(field, {'json_file': json_file, 'save_loc': save_loc}))

    return app
