function utc_time=now_utc()
% NOW_UTC works like "NOW", but returns the current time in UTC
% requires java, so won't work if you use matlab -nojava
utc_time = datenum(1970,1,1)+java.lang.System.currentTimeMillis/1000/3600/24;