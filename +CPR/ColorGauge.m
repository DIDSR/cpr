%% The 30-patch ColorGauge color target
% https://www.appliedimage.com/product-category/test-targets-and-charts/o-e-m-specific-test-charts/isa-colorgauge-targets/

classdef ColorGauge < CPR.ColorCheckerClass

    methods

        function obj = ColorGauge ()
            % define the color target

            obj.name = 'ColorGauge';
            obj.target_type = 2;

            obj.patch_n = 30;
            obj.all_patches = 1:obj.patch_n;
            obj.patch_row = 5;
            obj.patch_col = 6;
            obj.grayscale_patches = [8:11 14:17 20:23];                     % define the gray patches
            obj.chromatic_patches = obj.all_patches(~ismember(obj.all_patches,obj.grayscale_patches));  % the remaining ones are chromatic
            obj.white_patch = 8;                                            % the white patch

            obj.colorchecker = [...
                120   82   68;
                201  144  127;
                84  121  153;
                94  107   64;
                125  129  174;
                92  188  172;
                0  135  167;
                242  241  234;
                222  222  220;
                207  208  204;
                187  188  185;
                230  124   46;
                195   88  149;
                168  169  168;
                147  149  147;
                133  134  132;
                114  115  114;
                47   93  167;
                246  198   33;
                95   96   95;
                72   74   74;
                57   57   57;
                53   53   53;
                201   85   96;
                184   57   57;
                76  148   74;
                0   67  144;
                237  160   41;
                168  187   68;
                91   59  101];
        end

    end
end
