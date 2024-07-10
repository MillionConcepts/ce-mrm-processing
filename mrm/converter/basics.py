import datetime as dt


def modification_date():
    now = dt.datetime.now()
    return f"{now.year}-{str(now.month).zfill(2)}-{str(now.day).zfill(2)}"


def product_creation_time():
    return dt.datetime.now().isoformat()[:-7] + "Z"
