process hello {

    output:
    stdout

    script:
    """
    #!/usr/bin/python
    print("Hello World")
    """
}

workflow {
    // Run the hello process
    hello().view()
}

