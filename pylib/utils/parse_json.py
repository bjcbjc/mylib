import collections
import json
import re

__author__ = 'bjchen'


class NYGCTaskJson(object):
    def __init__(self, fn):
        self.process = json.load(open(fn))
        ## structure: [ pipeline, pipeline ...]
        ## pipeline is a dict with keys [pipeline_id, pipeline_name, tasks, pipeline_inputs, pipeline_outputs
        ## tasks is a list [task, task, ...]
        ## each task is a dict with keys [task_id, task_inputs, task_outputs, command_line, task_name]

    def get_pipeline_info(self, name, field):
        for pipeline in self.process:
            if pipeline['pipeline_name'] == name:
                for key in pipeline:
                    if key == 'pipeline_%s'%field:
                        return pipeline[key]
        return None

    def get_task_version(self):
        res = collections.OrderedDict()
        for pipeline in self.process:

            for task in pipeline['tasks']:
                name = task['task_name']
                version = task['version']
                if version is None: continue
                if name in res:
                    if res[name] != version:
                        raise ValueError('Multiple versions are used for the same task: ', name, version)
                else:
                    res[name] = version
        return res

    def get_command(self, abstract=False):
        name = list()
        command = list()
        for pipeline in self.process:
            for task in pipeline['tasks']:
                name.append(task['task_name'])
                cmd = str(task['command_line'])
                if abstract:
                    for i in xrange(len(task['task_inputs'])):
                        absFileName = self.abstract_file_name(task['task_inputs'][i])
                        cmd = cmd.replace(task['task_inputs'][i], '{%s}'%absFileName)
                        # cmd = cmd.replace(task['task_inputs'][i], '<input%02d>'%(i+1))
                    for i in xrange(len(task['task_outputs'])):
                        absFileName = self.abstract_file_name(task['task_outputs'][i])
                        cmd = cmd.replace(task['task_outputs'][i], '{%s}'%absFileName)
                        # cmd = cmd.replace(task['task_outputs'][i], '<output%02d>'%(i+1))
                    cmd = re.sub('-Djava.io.tmpdir=\S+', '-Djava.io.tmpdir={tmpdir}', cmd)
                    command.append(cmd)
                else:
                    command.append(task['command_line'])
        return name, command


    @staticmethod
    def abstract_file_name(fullPathFileName):
        baseName = fullPathFileName.split('/')[-1]
        sample = NYGCTaskJson.parse_sample_name(fullPathFileName)
        # determine if there is lane info in the file name
        laneInfo = re.findall('/Sample_%s/(%s\S+?)/'%(sample, sample), fullPathFileName)
        if len(laneInfo) > 0:
            laneInfo = laneInfo[0]
        if len(laneInfo) > 0 and laneInfo in baseName:
            res = baseName.replace(laneInfo, '').lstrip('_.')
        else:
            res = baseName.replace(sample, '').lstrip('_.')
        return res

    @staticmethod
    def parse_sample_name(fullPathFileName):
        tokens = list(set(re.findall('/Sample_(\S+?)/', fullPathFileName)))
        if len(tokens ) > 1:
            raise ValueError('more than one sample token found:', fullPathFileName, tokens)
        elif len(tokens) == 0:
            return ''
        return tokens[0]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-json', type=str, required=True, help='NYGC pipeline json')
    parser.add_argument('-action', type=str, choices=['get_exact_command', 'get_abstract_command', 'get_version', 'get_pipe_info'], required=True, help='command')
    parser.add_argument('-task', type=str, default=None, help='task name; if specified, only extract information related to the task')
    parser.add_argument('-pipeName', type=str, default=None, help='pipeline name; if specified, only extract <pipeField> information from the pipeline')
    parser.add_argument('-pipeField', type=str, default=None, help='pipeline field; eg. specify "inputs" to get pipeline_inputs under <pipeName>')
    parser.add_argument('-output', type=str, default=None, help='output file; stdout if not provided')

    args = vars(parser.parse_args())

    nygcJson = NYGCTaskJson(args['json'])

    if args['output'] is None:
        import sys
        out = sys.stdout
    else:
        out = open(args['output'], 'w')

    if args['action'] == 'get_exact_command' or args['action'] == 'get_abstract_command':
        name, command = nygcJson.get_command('abstract' in args['action'])
        for n, c in zip(name, command):
            if args['task'] is None or n == args['task']:
                out.write('##%s\n%s\n\n'%(n, c))
    elif args['action'] == 'get_version':
        res = nygcJson.get_task_version()
        for task, version in res.iteritems():
            if args['task'] is None or task == args['task']:
                out.write('%s\t%s\n'%(task, version))
    elif args['action'] == 'get_pipe_info':
        res = nygcJson.get_pipeline_info(args['pipeName'], args['pipeField'])
        if res is not None:
            for val in res:
                out.write('%s\n'%val)
    out.close()

