%% The 20-patch SFRPlus color target
% https://www.imatest.com/docs/sfrplus_instructions/
% Not recommended because it does not include the grayscale ramp

classdef SFRPlus < CPR.ColorCheckerClass

    methods

        function obj = SFRPlus ()
            % define the color target

            obj.name = 'SFRPlus';
            obj.target_type = 4;
            
            obj.patch_n = 20;
            obj.patch_row = 4;
            obj.patch_col = 5;
            obj.all_patches = [1:20];
            obj.grayscale_patches = [];                                     % there is no gray patch to evalute the gray ramp
            obj.chromatic_patches = obj.all_patches(~ismember(obj.all_patches,obj.grayscale_patches)); % all chromatic patch
            obj.white_patch = -1;                                           % there is no white patch

            obj.colorchecker = [...
                88	117	154;
                180	56	58;
                69	142	72;
                43	68	146;
                85	102	63;
                126	122	171;
                0	130	166;
                188	79	145;
                233	191	23;
                94	185	172;
                111	78	66;
                193	143	127;
                158	182	62;
                224	154	39;
                216	118	45;
                135	139	157;
                151	136	127;
                64	87	164;
                194	80	95;
                86	56	103];
        end

    end
end
