function [paths] = getobs2paths(getobsxml, outputdirectory)
    import com.asascience.ioos.netcdf.CreateNetcdf;
    import com.asascience.ioos.parser.GetObservationParser
    %import com.asascience.ioos.parser.SweDataRecordParser;
    %import ucar.nc2.ft.FeatureDatasetFactoryManager;
    %import ucar.nc2.dataset.NetcdfDataset;
    
    gop = com.asascience.ioos.parser.GetObservationParser();
    getObsModel = gop.parseGO(getobsxml);

    %memObs = getObsModel.getMemberObservation();
    %sweParser = SweDataRecordParser('timeSeriesProfile');
    %sweParser.parseSweDataRecord(singlestation, memObs.get(0));
     
    ncCreate  = CreateNetcdf(getObsModel);

    ncList = ncCreate.generateNetcdf(outputdirectory);
    for i = 1:length(ncList)
        filename = char(ncList(index));
        filename = strsplit(filename, '{');
        paths(i) = filename{1}(9:end-1);
        ncList.get(i-1).close()
    end
end
