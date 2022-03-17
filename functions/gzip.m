function files = gzip(my_directory)
    files = gunzip([my_directory,'\*.gz']);
end