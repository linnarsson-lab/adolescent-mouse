from typing import *
import cytograph as cg
import logging as lg
import luigi
import random
import string

lg.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=lg.DEBUG)


def logging(task: luigi.Task, log_dependencies: bool = False) -> lg.Logger:
	logger_name = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
	log_file = task.output().path + ".log"
	logger = lg.getLogger(logger_name)
	formatter = lg.Formatter('%(asctime)s %(levelname)s: %(message)s')
	fileHandler = lg.FileHandler(log_file, mode='w')
	fileHandler.setFormatter(formatter)
	# streamHandler = lg.StreamHandler()
	# streamHandler.setFormatter(formatter)

	logger.setLevel(lg.INFO)
	logger.addHandler(fileHandler)
	# logger.addHandler(streamHandler)

	if log_dependencies:
		logger.info("digraph G {")
		graph: Dict[str, List[str]] = {}

		def compute_task_graph(task: luigi.Task) -> None:
			name = task.__str__().split('(')[0]
			for dep in task.deps():
				if name in graph:
					graph[name].append(dep.__str__().split('(')[0])
				else:
					graph[name] = [dep.__str__().split('(')[0]]
				compute_task_graph(dep)

		compute_task_graph(task)
		for k, v in graph.items():
			for u in v:
				logger.info('"' + u + '" -> "' + k + '";')
		logger.info("}")
		logger.info("")

	for p in task.get_param_names():
		logger.info(f"{p} = {task.__dict__[p]}")
	logger.info("===")
	return logger
