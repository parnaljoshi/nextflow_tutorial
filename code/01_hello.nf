nextflow.enable.dsl=2

process hello {

    output:
    stdout

    script:
    """
    echo "Welcome to the world of Nextflow!"
    """
}

workflow {
    // Run the hello process
    hello().view()
}
