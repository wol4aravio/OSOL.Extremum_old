from argparse import ArgumentParser, RawTextHelpFormatter

import os
import shutil
import json

import numpy as np
import pandas as pd


def make_tasks(task_template, initials, u_min, u_max, location):
    start_config = []
    for p_init in np.linspace(initials['p']['min'], initials['p']['max'], initials['p']['N']).tolist():
        for q_init in np.linspace(initials['q']['min'], initials['q']['max'], initials['q']['N']).tolist():
            for r_init in np.linspace(initials['r']['min'], initials['r']['max'], initials['r']['N']).tolist():
                start_config.append([
                    {'name': 'p', 'value': p_init},
                    {'name': 'q', 'value': q_init},
                    {'name': 'r', 'value': r_init}
                ])
    start_config = dict([('{0:07d}'.format(id + 1), inits) for id, inits in enumerate(start_config)])

    if os.path.exists(location):
        shutil.rmtree(location)
    os.makedirs(location)

    for task_id, initial_conditions in start_config.items():
        temp_task = task_template.copy()
        temp_task['initial_conditions'] = initial_conditions
        for i in range(0, 3):
            temp_task['control_bounds'][i]['min'] = u_min
            temp_task['control_bounds'][i]['max'] = u_max
        json.dump(temp_task, open(location + '/{}.json'.format(task_id), 'w'), indent=4)

    legend = [[v['value'] for v in values] + [file_name] for file_name, values in start_config.items()]
    legend = pd.DataFrame(data=legend, columns=['p0', 'q0', 'r0', 'name'])
    legend.to_csv(location + '/legend.csv', index=False)
    return


initial_values = {
    'p': {'min': -25.0, 'max': 25.0, 'N': 51},
    'q': {'min': -25.0, 'max': 25.0, 'N': 51},
    'r': {'min': -25.0, 'max': 25.0, 'N': 51}
}


def main(args):
    make_tasks(task_template=json.load(open('TaskTemplateExplicit.json', 'r')),
               initials=initial_values, u_min=-500, u_max=500, location=args.folder_1)

    make_tasks(task_template=json.load(open('TaskTemplate.json', 'r')),
               initials=initial_values, u_min=-500, u_max=500, location=args.folder_2 )


parser = ArgumentParser(description='Task Creator', formatter_class=RawTextHelpFormatter)

parser.add_argument('--folder_1',
                    type=str,
                    help='Folder for explicit tasks')

parser.add_argument('--folder_2',
                    type=str,
                    help='Folder for final tasks')


args = parser.parse_args()

main(args)
