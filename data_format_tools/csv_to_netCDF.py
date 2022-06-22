
def csv_to_netCDF(infile,outfile, \
                    fill_value=-9999, \
                    outdimname='time' ):
    import netCDF4 as nc
    print( 'Converting '+infile+' to '+outfile)
    # open and read lines from input csv file
    inf=open(infile,'r')
    in_lines=inf.readlines()
    inf.close()
    
    # assume first line is headers
    headers = in_lines.pop(0)
    headers = headers[:-2]
    headers = headers.split(',')

    # Create Dictionary of lists to store data
    data_dict = { header:[] for header in headers }

    # loop over each data line storing value into correct place in dictionary
    for line in in_lines:
        split = line.split(',')
        for hdr,val in zip(headers,split):
            data_dict[hdr].append(float(val))

    # length of the dataset
    nt = len(data_dict[headers[0]])
    
    # open output file
    outf=nc.Dataset(outfile,'w')
    outf.createDimension(outdimname,nt)

    for hdr in headers:
        outvar=outf.createVariable(hdr,'float32',(outdimname),fill_value=fill_value)
        outvar[:]=data_dict[hdr]
    
    outf.close()


