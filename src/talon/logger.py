import logging


def _init_logger(verbosity):
    # https://coralogix.com/blog/python-logging-best-practices-tips/
    # https://stackoverflow.com/questions/14097061/easier-way-to-enable-verbose-logging

    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    level = levels[min(verbosity, len(levels) - 1)]  # cap to last level index

    # set defaults
    msg_fmt = "%(asctime)s : %(levelname)s : [%(filename)s:%(lineno)d] : %(message)s"
    date_fmt = "[ %Y-%m-%d %H:%M:%S ]"

    logging.basicConfig(level=level, format=msg_fmt, datefmt=date_fmt)
