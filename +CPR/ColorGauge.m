%% ColorChecker in 2D and AR/VR
classdef ColorGauge < CPR.ColorCheckerClass
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
    end

    methods

        function obj = ColorGauge ()

            obj.target_type = 2;
            obj.patch_n = 30;
            obj.all_patches = 1:obj.patch_n;
            obj.patch_row = 5;
            obj.patch_col = 6;
            obj.grayscale_patches = [8:11 14:17 20:23];
            obj.chromatic_patches = obj.all_patches(~ismember(obj.all_patches,obj.grayscale_patches));
            obj.white_patch = 8;

            %% probably RezChecker
            %             obj.patch_row = 6;
            %             obj.patch_col = 5;
            %             obj.grayscale_patches = [10:13 16:19 22:25];
            %             obj.chromatic_patches = [1:4 5 6 7 8 9 14 15 20 21 26 27:30];
            %             obj.white_patch = 8;

            % from K220154
            % obj.colorchecker = ones(30,3);

            obj.colorchecker = [...
                120.1090   82.4273   68.1429;
                201.6719  144.7946  127.1166;
                84.5178  121.7666  153.7429;
                94.9205  107.9609   64.5855;
                125.4773  129.9344  174.9751;
                92.0534  188.2957  172.0644;
                0  135.8458  167.4592;
                242.3976  241.0517  234.0347;
                222.7182  222.8318  220.2465;
                207.8550  208.3222  204.3309;
                187.5733  188.5992  185.7131;
                230.3748  124.9460   46.3150;
                195.7765   88.2185  149.4790;
                168.6420  169.9103  168.1632;
                147.8235  149.3388  147.5852;
                133.1231  134.5519  132.9758;
                114.3459  115.5798  114.5916;
                47.0046   93.2510  167.8959;
                246.3543  198.3829   33.1606;
                95.6482   96.7464   95.9883;
                72.9622   74.0468   74.0700;
                57.5282   57.0721   57.7471;
                53.5002   53.2531   53.5681;
                201.8058   85.2017   96.8953;
                184.0620   57.5330   57.5666;
                76.2764  148.9084   74.6631;
                0   67.3991  144.3738;
                237.9304  160.9509   41.8190;
                168.1059  187.2581   68.2507;
                91.2983   59.5369  101.2241];

        end

    end
end
