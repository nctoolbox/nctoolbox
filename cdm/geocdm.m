classdef geocdm < cdm


    methods

        function obj = geocdm(url)
            obj = obj@cdm(url)
        end

        function v = variable(obj, variableName)
            v = geocdmvariable(obj, variableName);
        end
                

    end


end