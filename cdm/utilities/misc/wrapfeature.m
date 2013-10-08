function dataset = wrapfeature(path)
    if ischar(path)
        netcdffile = ucar.nc2.NetcdfFile.open(path);
        netcdfdataset = ucar.nc2.dataset.NetcdfDataset(netcdffile);
    else
        netcdfdataset = path;
    end
    feature_type = ucar.nc2.ft.FeatureDatasetFactoryManager.findFeatureType(netcdfdataset);
    dataset = ucar.nc2.ft.FeatureDatasetFactoryManager.wrap(feature_type, netcdfdataset, [], java.util.Formatter());
end
