function setup_ncdataset_java
addjars(fileparts(which(mfilename)));

% configure log4j. Log4j is included in Matlab's classpath
root = org.apache.log4j.Logger.getRootLogger();
root.addAppender(org.apache.log4j.ConsoleAppender(org.apache.log4j.PatternLayout('%d{ISO8601} [%t] %-5p %c %x - %m%n')));
root.setLevel(org.apache.log4j.Level.WARN)