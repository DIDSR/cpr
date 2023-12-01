%% The 42-patch RezChenger color target
% https://www.imagescienceassociates.com/colorgauge-rezchecker-nano-target.html

classdef RezChecker < CPR.ColorCheckerClass

    methods

        function obj = RezChecker ()
            % define the color target

            obj.name = 'RezChecker';
            obj.target_type = 3;

            obj.patch_n = 42;
            obj.patch_row = 7;
            obj.patch_col = 6;
            obj.all_patches = [2:5 7 12 13 18 19:24 25:30 31:36 38:41];     % need to exclude empty patches
            obj.grayscale_patches = [20:23 26:29 32:35];                    % define the gray patches   
            obj.chromatic_patches = obj.all_patches(~ismember(obj.all_patches,obj.grayscale_patches)); % the remaining ones are chromatic
            obj.white_patch = 20;                                           % the white patch

            obj.colorchecker = [...
                  0    0    0;
                130   89   70;
                222  153  138;
                 69  126  162;
                102  110   72;
                  0    0    0;
                0  130  177;
                0  0         0;
                0  0         0;
                0  0         0;
                0  0         0;
                129  130  181;
                211  103  151;
                0  0         0;
                0  0         0;
                0  0         0;
                0  0      0;
                68  187  185;
                255  208   44;
                255  255  255;
                230  230  230;
                214  216  216;
                190  193  195;
                251  141   44;
                206   69   52;
                170  173  174;
                145  148  149;
                135  138  138;
                117  119  118;
                0    99  174;
                81  146   84;
                97   99   97;
                75   75   73;
                60   56   52;
                57   53   47;
                228  100  97;
                0    0    0;
                0    72  150;
                255  174   43;
                182  188   84;
                101   66  106;
                0    0    0];

        end

    end
end
