# ECP tools for manipulating datetime objects
import datetime
 
# create a incremental list of datetime objects
# requires a datetime.timedelta for tdelta
def DTarange(start_time,end_time,tdelta=datetime.timedelta(days=1)):
    datetime_list = [start_time]

    while datetime_list[-1]<end_time:
        datetime_list.append( datetime_list[-1]+tdelta )

    return datetime_list

# tool to add on an incremtn of months.
def add_months(sourcedate,months):
    import calendar
    month = sourcedate.month - 1 + months
    year = int(sourcedate.year + month / 12 )
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return datetime.datetime(year,month,day)

# create incremental list of datetime objects 
#  where increments are integer number of months
def DTarange_months(start_time,end_time,month_delta=1):
    datetime_list = [start_time]
    while datetime_list[-1]<end_time:
        datetime_list.append( add_months(datetime_list[-1],month_delta) )
    return datetime_list

