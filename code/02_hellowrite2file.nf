process hello {
    output:
    path 'result.txt'

    script:
        """
        echo "Welcome to the world of Nextflow!" > result.txt
        """
}

workflow {
    hello()
}
