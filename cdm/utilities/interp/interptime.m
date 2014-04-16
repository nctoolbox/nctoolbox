        
        function it = interptime(data, time, intime)
          % NCGEOVARIABLE.INTERPTIME - Interpolate data matrix to a particular time or time vector.
          % Useage: >> data = interptime(data, datatime, [2005 1 1 0 0 0]);
          %
          % NCTOOLBOX (https://github.com/nctoolbox/nctoolbox)
          ts = timeseries(data, time); % Create matlab timeseries
          wanttime = datenum(intime);
          tsnew = resample(ts, wanttime); % Resample to input time
          it = tsnew.data;
          
        end