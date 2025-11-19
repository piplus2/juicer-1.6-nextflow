process RUN_JAVA_TOOL {
    label "juicertools"

    debug true

    script:
    """
    export LC_ALL=en_US.UTF-8
    export _JAVA_OPTIONS="-Xmx${params.java_mem}"

    java -XshowSettings:vm -version
    java -XX:+PrintFlagsFinal -version | grep -i MaxHeapSize
    """
}
