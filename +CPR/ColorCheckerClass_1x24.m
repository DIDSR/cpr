%% ColorChecker in a linear shape

classdef ColorCheckerClass_1x24 < CPR.ColorCheckerClass
    %ColorCheckerClass_1x24 Rearrange ColorChecker as 1x24 

    methods

        function obj = ColorCheckerClass_1x24 
            % Constructor: define the layout
            
            obj.patch_col = 24;
            obj.patch_row = 1;
        end

    end
    
end