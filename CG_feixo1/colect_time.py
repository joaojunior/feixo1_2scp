import resource

def cpu_time():
    return resource.getrusage(resource.RUSAGE_SELF)[0]