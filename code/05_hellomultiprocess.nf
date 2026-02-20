nextflow.enable.dsl=2

process WRITE_TO_FILE {

	output:
	path "hello.txt"

	script:
	"""
	echo "Hello World from Process 1" >hello.txt
	"""
}

process PYTHON_HELLO {
	output:
	stdout

	script:
	"""
	#!/usr/bin/pytho
	print("Hello World from Process 2 (Python)")
	"""
}

process BASH_HELLO {

	output:
	stdout

	script:
	"""
	echo "Hello World from Process 3 (Bash)"
	"""
}

workflow {
	WRITE_TO_FILE()
	PYTHON_HELLO().view()
	BASH_HELLO().view()
}
