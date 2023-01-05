classdef ColorPerformanceReview < handle
    %COLORPERFORMANCEREVIEW Class for the Color Performance Review (CPR) tool
    %

    properties

        color_target          % color target
        all_patches           % list of all patches that are not empty 
        grayscale_patches     % list of grayscale patches
        chromatic_patches     % list of chromatic patches 
        all_patches_mask      % mask of all patches that are not empty 
        grayscale_patches_mask % mask of grayscale patches
        chromatic_patches_mask % mask of chromatic patches 

        devDataOriginal       % measurement data of the original (reference)
        devDataReproduced     % measurement data of the device (subject device)

        labOriginal       % LAB values of the original
        labReproduced     % LAB values of the device

        rgbOriginal       % RGB values of the original
        rgbReproduced     % RGB values of the device

        lchOriginal       % LCH values of the original
        lchReproduced     % LCH values of the device

        ilablchJoint      % joint table of id, L*, a*, b*, L*, c, h

        order_l           % order in lightness
        order_c           % order in chroma
        order_h           % order in hue
        order_a           % order in a*
        order_b           % order in b*

        im12_l6           % 1D order plot of the gray patches in lightness
        im12_l            % 1D order plot of the chromatic patches in lightness
        im12_c            % 1D order plot of the chromatic patches in chroma
        im12_h            % 1D order plot of the chromatic patches in hue

        kendall_l6        % Kendall tau-A of the gray patches in lightness
        kendall_l         % Kendall tau-A of the chromatic patches in lightness
        kendall_c         % Kendall tau-A of the chromatic patches in chroma
        kendall_h         % Kendall tau-A of the chromatic patches in hue

        target_type = 1;  % 24-patch, 30-patch
        patch_n           % number of patches, getting from CIELAB data

        CCE_name = [00 94 76];
        
        dE_from_ref            % dE00, dE94, and dEab
        dEmean                 % stats for dE00
        dEstd                  % stats for dE00
        dEmax                  % stats for dE00 
        dEmin                  % stats for dE00
        dEmedian               % stats for dE00  

        CCEVectorOriginal      % all dE from truth as a vector
        CCEVectorReproduced    % all dE from device as a vector

        CCE00TableOriginal     % all dE00 from truth as a table for the original
        CCE94TableOriginal     % all dE94 from truth as a table for the original
        CCEabTableOriginal     % all dEab from truth as a table for the original

        CCE00TableReproduced   % all dE00 from truth as a table for the reproduced
        CCE94TableReproduced   % all dE94 from truth as a table for the reproduced
        CCEabTableReproduced   % all dEab from truth as a table for the reproduced

        CCE_ratio              % CCE ratios for dE00, dE94, and dEab
        CCE_mean               % mean CCE ratios for dE00, dE94, and dEab

        datapath               % path for retrieving sample endoscopic images

        K_number               % K number     
        MakeName               % manufacturer name   
        DeviceName             % device name  
    end

    %% Constructor
    methods

        function obj = ColorPerformanceReview (devDataRef, devDataSub)

            [folder, name, ext] = fileparts(which('CPR.ColorPerformanceReview'));                 % where the data file is
            obj.datapath = [folder];                     % store the datapath

            switch size(devDataRef.Lab,1)
                case 24
                    obj.color_target = CPR.ColorCheckerClass;               % default ColorChecker
                case 30
                    obj.color_target = CPR.ColorGauge;               % default ColorChecker
                case 42
                    obj.color_target = CPR.RezChecker;               % default ColorChecker
                otherwise
                    disp '!!!Color Target Unknown!!!'
                    return
            end                        

            patchid = [1:obj.color_target.patch_n]';

            obj.all_patches = obj.color_target.all_patches';
            obj.grayscale_patches = obj.color_target.grayscale_patches';
            obj.chromatic_patches = obj.color_target.chromatic_patches';
            
            obj.all_patches_mask = ismember(patchid,obj.all_patches);
            obj.grayscale_patches_mask = ismember(patchid,obj.grayscale_patches);
            obj.chromatic_patches_mask = ismember(patchid,obj.chromatic_patches);

            assert(size(devDataRef.Lab,1)==size(devDataSub.Lab,1)); % check whether Ref and Sub have the same size

            obj.patch_n = size(devDataRef.Lab,1);                   % save the patch count

            assert(obj.patch_n == obj.color_target.patch_n);        % check whether Ref and ColorChecker have the same patch count  

            % save the original DeviceData
            obj.devDataOriginal = devDataRef;
            obj.devDataReproduced = devDataSub;

            % save the LAB
            obj.labOriginal = devDataRef.Lab;
            obj.labReproduced = devDataSub.Lab;
            
            % calculate the RGB
            obj.rgbOriginal = lab2rgb(obj.labOriginal)*255;
            obj.rgbReproduced = lab2rgb(obj.labReproduced)*255;

            obj.rgbOriginal = min(obj.rgbOriginal,255);
            obj.rgbOriginal = max(obj.rgbOriginal,0);
            
            obj.rgbReproduced = min(obj.rgbReproduced,255);
            obj.rgbReproduced = max(obj.rgbReproduced,0);

            % calculate the LCH
            obj.lchOriginal = obj.lab2lch(obj.labOriginal);
            obj.lchReproduced = obj.lab2lch(obj.labReproduced);

            obj.sort_labch;

            obj.calcualte_concordance;

            obj.calculate_dE_from_truth;

            obj.calculate_CCE;

            obj.calculate_visual_1x24;
        end

        function sort_labch (obj)
            % Sort patches by L*, a*, b*, c, and h

            patch_n = obj.patch_n;

            % construct ilabch
            ilablchOriginal = [[1:patch_n]' obj.labOriginal obj.lchOriginal];
            ilablchReproduced = [[1:patch_n]' obj.labReproduced obj.lchReproduced];

            ilablchJoint = [ilablchOriginal ilablchReproduced];

            obj.ilablchJoint = ilablchJoint;

            obj.order_l = helper(ilablchJoint,2,9);
            obj.order_a = helper(ilablchJoint,3,10);
            obj.order_b = helper(ilablchJoint,4,11);
            obj.order_c = helper(ilablchJoint,6,13);
            obj.order_h = helper(ilablchJoint,7,14);

            return

            function order = helper (tab,col_ref,col_subject)
                % Sort "tab" by col_ref and by col_subject

                sorted_for_ref = sortrows(tab,col_ref);
                order(:,1) = sorted_for_ref(:,1);

                sorted_for_subject = sortrows(tab,col_subject);
                order(:,2) = sorted_for_subject(:,8);
            end
        end

        function calcualte_concordance (obj)
            % Calculate kendall_l6, kendall_l, kendall_c, and kendall_h

            x = obj.labOriginal(obj.grayscale_patches_mask,1);
            y = obj.labReproduced(obj.grayscale_patches_mask,1);
            obj.kendall_l6 = obj.kendall_tau_a(x,y);

            x = obj.labOriginal(obj.chromatic_patches_mask,1);
            y = obj.labReproduced(obj.chromatic_patches_mask,1);
            obj.kendall_l = obj.kendall_tau_a(x,y);

            x = obj.lchOriginal(obj.chromatic_patches_mask,2);
            y = obj.lchReproduced(obj.chromatic_patches_mask,2);
            obj.kendall_c = obj.kendall_tau_a(x,y);

            x = obj.lchOriginal(obj.chromatic_patches_mask,3);
            y = obj.lchReproduced(obj.chromatic_patches_mask,3);
            obj.kendall_h = obj.kendall_tau_a(x,y);

        end

        function calculate_dE_from_truth (obj)

            dE = zeros(obj.patch_n,3);
            
            % mask = [1:obj.white_patch-1 obj.white_patch+1:obj.patch_n];  % remove #19 white??
            mask = obj.all_patches;

            for i = 1:obj.patch_n

                [dE00 dE94 dEab] = CPR.ColorClass.LAB2dE(obj.labOriginal(i,:)',obj.labReproduced(i,:)');
                
                dE(i,1:3) = [dE00 dE94 dEab];

            end

            obj.dE_from_ref = dE;

            obj.dEmean = mean(dE(mask));
            obj.dEstd = std(dE(mask));
            obj.dEmax = max(dE(mask));
            obj.dEmin = min(dE(mask));
            obj.dEmedian = median(dE(mask));
        end

        function calculate_CCE (obj)
            % Analyze CCE

            devOriginal = obj.devDataOriginal;
            devReproduced = obj.devDataReproduced;

            n = numel(obj.all_patches);

            truth(1:n) = CPR.ColorClass([0 0 0]);
            measured(1:n) = CPR.ColorClass([0 0 0]);

            for i = 1:n
                ii = obj.all_patches(i);
                truth(i) = devOriginal.Lab(ii,:);
                measured(i) = devReproduced.Lab(ii,:);
            end

            % dE pairwise
            obj.CCEVectorOriginal = zeros(n*(n-1)/2,3);
            obj.CCEVectorReproduced = zeros(n*(n-1)/2,3);

            obj.CCE00TableOriginal = zeros(n,n);
            obj.CCE94TableOriginal = zeros(n,n);
            obj.CCEabTableOriginal = zeros(n,n);

            obj.CCE00TableReproduced = zeros(n,n);
            obj.CCE94TableReproduced = zeros(n,n);
            obj.CCEabTableReproduced = zeros(n,n);

            % the color patches considered for evaluation
            k = 0;
            for ii = 1:n
                for jj = ii+1:n
                    k = k + 1;
                    i = ii;
                    j = jj;

                    % do reference
                    [dE00 dE94 dEab] = CPR.ColorClass.LAB2dE (truth(i).lab', truth(j).lab');
                    obj.CCEVectorOriginal(k,:) = [dE00 dE94 dEab]';

                    obj.CCE00TableOriginal(i,j) = dE00;
                    obj.CCE94TableOriginal(i,j) = dE94;
                    obj.CCEabTableOriginal(i,j) = dEab;

                    % do device
                    [dE00 dE94 dEab] = CPR.ColorClass.LAB2dE (measured(i).lab', measured(j).lab');
                    obj.CCEVectorReproduced(k,:) = [dE00 dE94 dEab]';

                    obj.CCE00TableReproduced(i,j) = dE00;
                    obj.CCE94TableReproduced(i,j) = dE94;
                    obj.CCEabTableReproduced(i,j) = dEab;
                end
            end

            obj.CCE_ratio = obj.CCEVectorReproduced ./ obj.CCEVectorOriginal;
            obj.CCE_mean = mean(obj.CCE_ratio);
        end

        function calculate_visual_1x24 (obj)
            % Order in lightness, chroma, and hue - 1D view

            % 'Lightness6'
            obj.im12_l6 = helper(obj.order_l,obj.grayscale_patches,obj.kendall_l6);

            % 'Lightness'
            obj.im12_l = helper(obj.order_l,obj.chromatic_patches,obj.kendall_l);

            % 'Chroma'
            obj.im12_c = helper(obj.order_c,obj.chromatic_patches,obj.kendall_c);

            % 'Hue'
            obj.im12_h = helper(obj.order_h,obj.chromatic_patches,obj.kendall_h);

            return

            function im12 = helper (order,range,tau)
                cc1 = CPR.ColorCheckerClass_1x24;

                mask = ismember(order(:,1),range);
                idx1 = order(mask,1);
                rgb_l24 = obj.rgbOriginal(idx1,:);
                im1 = cc1.getImage_with_index(rgb_l24,idx1);

                mask = ismember(order(:,2),range);
                idx2 = order(mask,2);
                rgb_l24 = obj.rgbReproduced(idx2,:);
                im2 = cc1.getImage_with_index(rgb_l24,idx2);

                im12 = obj.stitching_two_colorbars(im1,im2,idx1,idx2);

                %
                % add tau
                %
                tau_str = sprintf('Tau %.2f',tau);
                im12 = insertText(im12,[10 size(im12,1)/2],tau_str,TextColor='black',FontSize=24,BoxColor='green');

            end

        end

        function order_1D_l6 (obj,ha)
            axes(ha)
            obj.layout_1D(obj.order_l,obj.grayscale_patches,obj.kendall_l6);
        end

        function order_1D_l (obj,ha)
            axes(ha)
            obj.layout_1D(obj.order_l,obj.chromatic_patches,obj.kendall_l);
        end

        function order_1D_c (obj,ha)

            % 'Chroma'
            axes(ha)
            obj.layout_1D(obj.order_c,obj.chromatic_patches,obj.kendall_c);
        end

        function order_1D_h (obj,ha)
            % 'Hue'
            axes(ha)
            obj.layout_1D(obj.order_h,obj.chromatic_patches,obj.kendall_h);
        end

        function layout_1D (obj,order,range,tau)
            cc1 = CPR.ColorCheckerClass_1x24;

            mask = ismember(order(:,1),range);
            idx1 = order(mask,1);
            rgb1 = obj.rgbOriginal(idx1,:);

            mask = ismember(order(:,2),range);
            idx2 = order(mask,2);
            rgb2 = obj.rgbReproduced(idx2,:);

            offset = 0.4;
            row_from = 0.5+offset;
            row_to = 0.5-offset;
            markersize = 15;
            col_gap = 1/20;

            %hf = gca;
            %hf.Units = 'Inches';

            hold on

            for i = 1:numel(idx1)
                mycolor = rgb1(i,:)/255;
                mystr = sprintf('%02d',idx1(i));

                xsize = 1/21;
                ysize = 1/5;
                pos = [i*col_gap row_from xsize ysize];

                %text(pos(1),pos(2),mystr);
                %upper{i} = annotation('rectangle',pos,'FaceColor',mycolor);
                plot(pos(1),pos(2),'s','MarkerSize',markersize,'MarkerFaceColor',mycolor,'MarkerEdgeColor',mycolor);
                upper{i} = pos;
            end

            for i = 1:numel(idx2)
                mycolor = rgb2(i,:)/255;
                mystr = sprintf('%02d',idx2(i));
                pos = [i*col_gap row_to xsize ysize];
%                 ha = annotation('textbox',pos, ...
%                     'String',mystr,...
%                     'BackgroundColor',mycolor);
                % lower{i} = annotation('rectangle',pos,'FaceColor',mycolor);
                plot(pos(1),pos(2),'s','MarkerSize',markersize,'MarkerFaceColor',mycolor,'MarkerEdgeColor',mycolor);
                lower{i} = pos;
            end

            for location1 = 1:numel(idx1)
                i_from = idx1(location1);
                i_to = idx2(location1);

                patch_id = idx1(location1);
                location2 = find(idx2 == patch_id);
                mycolor = rgb1(location1,:)/255;

%                 x1 = upper{location1}.Position(1) + upper{location1}.Position(3)/2 ;
%                 y1 = upper{location1}.Position(2);
%                 x2 = lower{location2}.Position(1) + upper{location1}.Position(3)/2;
%                 y2 = lower{location2}.Position(2) + ysize;
  
                x1 = upper{location1}(1);
                y1 = upper{location1}(2);
                x2 = lower{location2}(1);
                y2 = lower{location2}(2);
                plot([x1 x2],[y1 y2],'-','Color',mycolor,'LineWidth',2);

                % ha = annotation('line',[x1 x2],[y1 y2]);
                %ha.Parent = hf.CurrentAxes;
                
            end

            axis off
            axis([0 1 0 1.2])
            %im12 = obj.stitching_two_colorbars(im1,im2,idx1,idx2);

            %
            % add tau
            %
            %tau_str = sprintf('Tau %.2f',tau);
            %im12 = insertText(im12,[10 size(im12,1)/2],tau_str,TextColor='black',FontSize=24,BoxColor='green');

        end
        
    end

    %% Functions used by the class
    methods

        function tau_a = kendall_tau_a (obj,x,y)
            % Kendall rank correlation coefficient
            % https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient

            n = length(x);
            pairs = 0;
            concordant = 0;
            discordant = 0;
            for i = 1:n
                for j = i+1:n
                    pairs = pairs + 1;
                    if sign(x(i)-x(j)) * sign(y(i)-y(j)) >= 0
                        concordant = concordant+1;
                    else
                        discordant = discordant+1;
                    end
                end
            end
            %[concordant discordant pairs]
            tau_a = (concordant-discordant)/pairs;
        end

        function [p, yfit, rsq] = linear_regression (obj,x,y)
            % linear regression of a scatter plot

            p = polyfit(x,y,1);
            yfit = polyval(p,x);
            yresid = y - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(y)-1) * var(y);
            rsq = 1 - SSresid/SStotal;
        end

        function lch = lab2lch (obj,lab)
            % Convert CIELAB to CIELCH

            l = lab(:,1);
            a = lab(:,2);
            b = lab(:,3);
            lch(:,1) = l;
            lch(:,2) = (a.^2 + b.^2).^0.5;
            lch(:,3) = atan2d(b,a);
        end

        function save_figure (obj,resizeto,filename)
            % Save the plot with the given size

            my_gcf = gcf;
            my_gcf.Units = 'pixels';
            my_gcf.Position = resizeto;
            saveas(gcf,['fig/' filename])
        end
        
    end

    %% Functions using the class
    methods

        function evaluate_dE_from_truth (obj,ha1,ha2)
            metric = 1;

            axes(ha1)
            h1 = bar(obj.all_patches,obj.dE_from_ref(obj.all_patches_mask,metric),'FaceColor','flat');
            h1.CData(:,:) = obj.rgbOriginal(obj.all_patches_mask,:)/255;
            xlabel('Patch #')
            ylabel('Color Difference \DeltaE_{00}')
            title(sprintf('\\mu=%0.2f, \\delta=%0.2f', ...
                obj.dEmean(metric), ...
                obj.dEstd(metric)))

            axes(ha2)
            set(gca,'TickLabelInterpreter','latex');

%            boxplot(obj.dE_from_ref(obj.all_patches_mask,metric))
            boxplot(obj.dE_from_ref(obj.all_patches_mask,1:3),{'dE00','dE94','dE76'});

            ylabel('Color Difference \DeltaE_{00}')
            axis([0 4 h1.Parent.YLim(1) h1.Parent.YLim(2)])
            title(sprintf('min=%0.2f, median=%0.2f, max=%0.2f', ...
                obj.dEmin(metric), ...
                obj.dEmedian(metric), ...
                obj.dEmax(metric)))
        end

        function evaluate_CCE (obj,metric,ha)
            % Show dE_in vs. dE_out
            % metric: 1=dE00, 2=dE94, 3=dEab

            axes(ha)
            
            x = obj.CCEVectorOriginal(:,metric);
            y = obj.CCEVectorReproduced(:,metric);
            [p,yfit,rsq] = obj.linear_regression(x,y);
        
            yoverx = y./x;
            cce_good = nnz(yoverx >= 1);
            cce_all = numel(yoverx);
            cce_good_percent = cce_good / cce_all * 100;

            box = [0 100 0 100];

            hold on

            switch metric
                case 1
                    CCETableIn = obj.CCE00TableOriginal;
                    CCETableOut = obj.CCE00TableReproduced;
                case 2
                    CCETableIn = obj.CCE94TableOriginal;
                    CCETableOut = obj.CCE94TableReproduced;
                case 3
                    CCETableIn = obj.CCEabTableOriginal;
                    CCETableOut = obj.CCEabTableReproduced;
            end

            n = numel(obj.all_patches);
            for i = 1:n
                for j = i+1:n
                    dein = CCETableIn(i,j);
                    deout = CCETableOut(i,j);
                    colorin = obj.rgbOriginal(i,:)/255;
                    colorout = obj.rgbReproduced(i,:)/255;
                    plot(dein,deout,'_','MarkerEdgeColor',colorin);
                    plot(dein,deout,'|','MarkerEdgeColor',colorout);
                end
            end

            % plot(x,y,'o')

            plot(box(1:2),box(3:4),'--k')   % add the identity line

            %
            % add fitting data
            %
            if 0
                yfit = polyval(p,box(1):box(2));
                plot(box(1):box(2),yfit,':r')   % add the fitting line

                fitting_eq = sprintf('y=%.2fx+%.2f',p(1),p(2));
                text(box(1)+(box(2)-box(1))/10,box(3)+(box(4)-box(3))/10*9,fitting_eq);
            end

            %
            % add good CCE percent
            %
            if 1
                cce_eq = sprintf('%.2f%%',cce_good_percent);
                text(box(1)+(box(2)-box(1))/10,box(3)+(box(4)-box(3))/10*9,cce_eq);
            end

            grid on
            axis square
            axis(box)
            xlabel('Input \Delta{E}')
            ylabel('Output \Delta{E}')
            title(sprintf('CCE_{%02d}=%.2f%c%0.2f', ...
                obj.CCE_name(metric), ...
                obj.CCE_mean(metric), ...
                char(177),...
                std(obj.CCE_ratio(:,metric))))
            %legend('Measured','Identity')
            %legend('Location','Southeast')
        end

        function [im1, im2] = evaluate_visual_XYZ_im (obj)
            % Visualize CIEXYZ data

            XYZ_d65 = rgb2xyz([1 1 1]);
            Y_d65 = XYZ_d65(2);

            XYZ1_orig = obj.devDataOriginal.XYZ;
            Y1_orig = XYZ1_orig(19,2);
            XYZ1 = XYZ1_orig / Y1_orig * Y_d65;

            XYZ2_orig = obj.devDataReproduced.XYZ;
            Y2_orig = XYZ2_orig(19,2);
            XYZ2 = XYZ2_orig / Y2_orig * Y_d65;

            %
            %
            %
            ccc = CPR.ColorCheckerClass;

            rgb1 = xyz2rgb(XYZ1)*255;
            im1 = ccc.getImage(rgb1);

            rgb2 = xyz2rgb(XYZ2)*255;
            im2 = ccc.getImage(rgb2);
        end
        
        function evaluate_visual_XYZ (obj)
            [im1, im2] = evaluate_visual_XYZ_im (obj);

            subplot(1,2,1)
            imshow(im1)
            title('Reference')

            subplot(1,2,2)
            imshow(im2)
            title('Reproduced')
        end

        function [im0, im1, im2] = evaluate_visual_im (obj)
            ccc = obj.color_target;

            im0 = ccc.getImage_with_index;
            im1 = ccc.getImage(obj.rgbOriginal);
            im2 = ccc.getImage(obj.rgbReproduced);
        end

        function evaluate_visual (obj, ha0, ha1, ha2)
            % Visualize CIEXYZ data

            [im0, im1, im2] = obj.evaluate_visual_im;
            
            axes(ha0)
            imshow(im0)
            title('Original ColorChecker')

            axes(ha1)
            imshow(im1)
            title('Reference Data, LAB')

            axes(ha2)
            imshow(im2)
            title('Reproduced Data, LAB')
        end

        function im2 = simulate_im (obj,im1)
            % Interpolate the image to be simulated

            rgb1 = reshape(im1,size(im1,1)*size(im1,2),3);
            lab1 = rgb2lab(rgb1);

            labref = obj.labOriginal;
            labtarget = obj.labReproduced;

            Fl = scatteredInterpolant(labref,labtarget(:,1));
            Fa = scatteredInterpolant(labref,labtarget(:,2));
            Fb = scatteredInterpolant(labref,labtarget(:,3));

            mymethod = 'linear';

            Fl.Method = mymethod;
            Fa.Method = mymethod;
            Fb.Method = mymethod;

            lab2 = [Fl(lab1(:,1),lab1(:,2),lab1(:,3)) Fa(lab1(:,1),lab1(:,2),lab1(:,3)) Fb(lab1(:,1),lab1(:,2),lab1(:,3))];
            rgb2 = lab2rgb(lab2);
            im2 = reshape(rgb2,size(im1,1),size(im1,2),3);
        end

        function simulate_axes (obj,im,ha)
            axes(ha);
            imshow(im);
        end

        function simulate (obj,filename)
                        
            im1 = imread(filename);
            ha1 = subplot(1,2,1);

            im2 = obj.simulate_im (im1);
            ha2 = subplot(1,2,2);

            obj.simulate_axes(im1,ha1);
            title('Ground Truth')

            obj.simulate_axes(im2,ha2);
            title('Reproduced')

        end

        function im12 = stitching_two_colorbars (obj,im1,im2,idx1,idx2)
            % Stitch two colorbars together with lines

            im3 = im1;

            im3 = im3*0.01;

            im12 = [im1 ; im3; im2];

            n = length(idx1);

            idx2two = [idx2 [1:n]'];
            idx2twosorted = sortrows(idx2two,1);

            imw = size(im12,2);
            imh = size(im12,1);

            boxw = imw/n;

            posx = [1:length(idx1)]*boxw - boxw/2 ;
            posy = [boxw*1 boxw*2];

            for location1 = 1:n
                patch_id = idx1(location1);
                location2 = find(idx2 == patch_id);

                im12 = insertShape(im12,'Line',[posx(location1) posy(1) posx(location2) posy(2)],'Color','Yellow','LineWidth',5);
            end

            %             clf
            %             imshow(im12)
        end

        function order = evaluate_order_lightness_6 (obj, ha)
            % show linearity plot for lightness 6 patches

            axes(ha);

            l_original = obj.labOriginal(obj.grayscale_patches_mask,1);
            l_reproduced = obj.labReproduced(obj.grayscale_patches_mask,1);
            id = obj.grayscale_patches;

            [p, yfit, rsq] = obj.linear_regression(l_original,l_reproduced);

            obj.linearity_plot(l_original,l_reproduced,id,[0 100 0 100],p);

            xlabel('L* reference')
            ylabel('L* measured')
            title('Lightness Order')

        end
        
        function order = evaluate_order_lightness (obj, ha)
            % show linearity plot for lightness
    
            axes(ha);

            l_original = obj.labOriginal(obj.chromatic_patches_mask,1);
            l_reproduced = obj.labReproduced(obj.chromatic_patches_mask,1);
            id = obj.chromatic_patches;

            [p, yfit, rsq] = obj.linear_regression(l_original,l_reproduced);

            obj.linearity_plot(l_original,l_reproduced,id,[0 100 0 100],p);

            xlabel('L* reference')
            ylabel('L* measured')
            title('Lightness Order')

        end

        function linearity_plot (obj, x, y, id, box, p)
            % show one linearity plot

            hold on
            for i=1:length(id)
                rgb256 = obj.rgbOriginal(id(i),:);
                rgb = rgb256 / 255.0;
                rgb = min(1,rgb);
                rgb = max(0,rgb);
                plot(x(i),y(i),'o','MarkerFaceColor',rgb,'MarkerEdgeColor',rgb)
            end
            text(x,y+1,num2str(id))

            axis(box)
            plot(box(1:2),box(3:4),'-k')   % add the identity line
            yfit = polyval(p,box(1):box(2));
            plot(box(1):box(2),yfit,':r')   % add the fitting line

            fitting_eq = sprintf('y=%.2fx+%.2f',p(1),p(2));
            text(box(1)+(box(2)-box(1))/10,box(3)+(box(4)-box(3))/10*9,fitting_eq)

            axis square
            grid on            
        end

        function order = evaluate_order_chroma (obj,ha)
            % show linearity plot for chroma

            axes(ha);

            c_original = obj.lchOriginal(obj.chromatic_patches_mask,2);
            c_reproduced = obj.lchReproduced(obj.chromatic_patches_mask,2);
            id = obj.chromatic_patches;

            [p, yfit, rsq] = obj.linear_regression(c_original,c_reproduced);

            obj.kendall_c = obj.kendall_tau_a(c_original,c_reproduced);

            obj.linearity_plot(c_original,c_reproduced,id,[0 100 0 100],p);

            xlabel('Chroma reference')
            ylabel('Chroma measured')
            title('Chroma Order')

        end

        function order = evaluate_order_hue (obj,ha)
            % show linearity plot for hue

            axes(ha);

            h_original = obj.lchOriginal(obj.chromatic_patches_mask,3);
            h_reproduced = obj.lchReproduced(obj.chromatic_patches_mask,3);
            id = obj.chromatic_patches;

            [p, yfit, rsq] = obj.linear_regression(h_original,h_reproduced);

            obj.kendall_h = obj.kendall_tau_a(h_original,h_reproduced);

            obj.linearity_plot(h_original,h_reproduced,id,[-180 180 -180 180],p);
            
            xlabel('Hue reference')
            ylabel('Hue measured')
            title('Hue Order')

        end

        function show2dquiver_polar_select (obj,mask,ha)
            % Show hue/chroma transfer of a subset in CIELAB on a polar

            axes(ha)

            show3dquiver_select (obj,mask,ha)    
            
            line_color = 0.8;

            plot_circle(25)
            plot_circle(50)
            plot_circle(75)
            plot_circle(100)

            for i = 0:11
                th = 2*pi/12 * i;
                plot_radial(th,100)
            end

            grid off

            axs = axis;
            x1 = axs(1);
            x2 = axs(2);
            y1 = axs(3);
            y2 = axs(4);
            plot([x1 x2],[(y1+y2)/2 (y1+y2)/2],'-','Color',[1 1 1]*line_color);
            plot([(x1+x2)/2 (x1+x2)/2],[y1 y2],'-','Color',[1 1 1]*line_color);

            return


            function plot_circle (r)
                th = 0:0.01:2*pi;
                plot(r*cos(th),r*sin(th),'-','Color',[1 1 1]*line_color)
            end

            function plot_radial (th,r)
                plot([0 r*cos(th)],[0 r*sin(th)],'-','Color',[1 1 1]*line_color)
            end
        end

        function show3dquiver_select (obj,mask,ha)
            % Show color transfer of a subset in CIELAB
            
            axes(ha)

            % set mask to all if the argument is missing
%             if nargin == 3
%                 mask = 1:size(obj.labOriginal,1);
%             end

            lab1 = obj.labOriginal(mask,:);
            lab2 = obj.labReproduced(mask,:);

            rgb1 = lab2rgb(lab1);
            rgb1 = min(1,rgb1);
            rgb1 = max(0,rgb1);

            rgb2 = lab2rgb(lab2);
            rgb2 = min(1,rgb2);
            rgb2 = max(0,rgb2);

            dif = lab2-lab1;

            hold on
            quiver3(lab1(:,2),lab1(:,3),lab1(:,1),dif(:,2),dif(:,3),dif(:,1),0,'k.')
            for i = 1:length(mask)
                plot3(lab1(i,2),lab1(i,3),lab1(i,1),'o','MarkerFaceColor',rgb1(i,:),'MarkerEdgeColor',rgb1(i,:))
                plot3(lab2(i,2),lab2(i,3),lab2(i,1),'x','MarkerFaceColor',rgb2(i,:),'MarkerEdgeColor',rgb2(i,:))
            end
            zlabel('L*')
            xlabel('a*')
            ylabel('b*')
            grid on
            axis equal
        end

    end

    %% Funcitons for CPR Live Script
    methods

        function one_pager (obj)
            cpr = obj;

            clf
            hg = gcf;
            set(hg,'Visible','on')

            td = tiledlayout(8,6,'TileSpacing','tight');
            title(td,sprintf('%s %s %s',obj.K_number,obj.MakeName,obj.DeviceName))

            ha0 = nexttile(1,[2 1]);
            if strcmp(obj.devDataOriginal.source,'XYZ')
                [im1, im2] = obj.evaluate_visual_XYZ_im;
                imshow(im2)
                title('Reproduced Data, XYZ')
            else
                axis off
            end

            ha0 = nexttile(2,[2 1]);
            ha1 = nexttile(13,[2 1]);
            ha2 = nexttile(14,[2 1]);

            cpr.evaluate_visual(ha0,ha1,ha2);
            filename = [obj.datapath '/sample_polyp.png'];

            im1 = imread(filename);
            im2 = cpr.simulate_im (im1);

            ha1 = nexttile(25,[2 1]);
            ha2 = nexttile(26,[2 1]);

            cpr.simulate_axes(im1,ha1);
            title('Original Scene')

            cpr.simulate_axes(im2,ha2);
            title('Predicted Device Output')

            ha1 = nexttile(37,[2 1]);
            ha2 = nexttile(38,[2 1]);
            cpr.evaluate_dE_from_truth(ha1,ha2);

            ha1 = nexttile(27,[1 2]); cpr.order_1D_l6(ha1); title(sprintf('Lightness order of gray patches,{\\tau}=%.2f',cpr.kendall_l6))
            ha2 = nexttile(33,[1 2]); cpr.order_1D_l(ha2); title(sprintf('Lightness order of chromatic patches,{\\tau}=%.2f',cpr.kendall_l))
            ha3 = nexttile(39,[1 2]); cpr.order_1D_h(ha3); title(sprintf('Hue order of chromatic patches,{\\tau}=%.2f',cpr.kendall_h))
            ha4 = nexttile(45,[1 2]); cpr.order_1D_c(ha4); title(sprintf('Chroma order of chromatic patches,{\\tau}=%.2f',cpr.kendall_c))

%            ha1 = nexttile(27,[1 2]); cpr.check_order_axes(cpr.im12_l6,ha1); title(sprintf('Lightness order of gray patches,{\\tau}=%.2f',cpr.kendall_l6))
%            ha2 = nexttile(33,[1 2]); cpr.check_order_axes(cpr.im12_l,ha2); title(sprintf('Lightness order of chromatic patches,{\\tau}=%.2f',cpr.kendall_l))
%            ha3 = nexttile(39,[1 2]); cpr.check_order_axes(cpr.im12_h,ha3); title(sprintf('Hue order of chromatic patches,{\\tau}=%.2f',cpr.kendall_h))
%            ha4 = nexttile(45,[1 2]); cpr.check_order_axes(cpr.im12_c,ha4); title(sprintf('Chroma order of chromatic patches,{\\tau}=%.2f',cpr.kendall_c))

            ha1 = nexttile(3,[2 1]);
            cpr.evaluate_order_lightness_6(ha1); title(sprintf('Lightness order of gray patches'))

            ha2 = nexttile(4,[2 1]);
            cpr.evaluate_order_lightness(ha2); title(sprintf('Lightness order of chromatic patches'))

            ha3 = nexttile(15,[2 1]);
            cpr.evaluate_order_hue(ha3); title(sprintf('Hue order of chromatic patches'))

            ha4 = nexttile(16,[2 1]);
            cpr.evaluate_order_chroma(ha4); title(sprintf('Chroma order of chromatic patches'))


            ha1 = nexttile(5,[2 1]); cpr.show3dquiver_select(cpr.chromatic_patches,ha1); view([12 32]); axis equal; title('Color transfer')
            ha2 = nexttile(6,[2 1]); cpr.show2dquiver_polar_select(cpr.chromatic_patches,ha2); view([0 90]); axis([-100 100 -100 100 0 100]); title('Hue Chroma')
            ha3 = nexttile(17,[2 1]); cpr.show3dquiver_select(cpr.chromatic_patches,ha3); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('L* a* plane')
            ha4 = nexttile(18,[2 1]); cpr.show3dquiver_select(cpr.chromatic_patches,ha4); view([90 0]); axis square; axis([-100 100 -100 100 0 100]); title('L* b* plane')
            ha5 = nexttile(29,[2 1]); cpr.show3dquiver_select(cpr.grayscale_patches,ha5); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('Lightness gray')
            ha6 = nexttile(30,[2 1]); cpr.show2dquiver_polar_select(cpr.grayscale_patches,ha6); view([0 90]); axis square; axis([-30 30 -30 30 0 100]); title('Chromaticity gray')

            ha1 = nexttile(41,[2 1]); cpr.evaluate_CCE(1,ha1);
            % ha2 = nexttile(24); cpr.evaluate_CCE(2,ha2);
            % ha3 = nexttile(27); cpr.evaluate_CCE(3,ha3);

            ha1 = nexttile(42,[2 1]);
            boxplot(cpr.CCE_ratio,{'CCE00','CCE94','CCE76'})
            ylabel('CCE (%)')
            title('CCE by different \DeltaE')

            set(hg,'Position',[137 162 1714 939])
            saveas(hg,'one_pager.png')

        end

        function check_Lab_data (obj)
            % Visualize CIEXYZ data

            Ref_L_star = obj.devDataOriginal.Lab(:,1);
            Ref_a_star = obj.devDataOriginal.Lab(:,2);
            Ref_b_star = obj.devDataOriginal.Lab(:,3);

            Sub_L_star = obj.devDataReproduced.Lab(:,1);
            Sub_a_star = obj.devDataReproduced.Lab(:,2);
            Sub_b_star = obj.devDataReproduced.Lab(:,3);

            CIELAB_Data = table(Ref_L_star,Ref_a_star,Ref_b_star,Sub_L_star,Sub_a_star,Sub_b_star)
        end

        function check_Lab_visual (obj)
            % Visualize CIELAB data
            % The following charts show the simulated visual results when using the provided reference white. Use these charts to assess how the device would reproduce the ColorChecker.

            clf
            h_visual = gcf;

            ha0 = subplot(1,3,1);
            ha1 = subplot(1,3,2);
            ha2 = subplot(1,3,3);            
            obj.evaluate_visual(ha0,ha1,ha2);
            set(h_visual,'Units','normalized','Position',[0 0 0.5 0.25]);
            obj.save_figure([48 267 676 229],'labx2.png')
        end

        function check_endoscopic_scene (obj,simulate_filename)
            % Visualize endoscopic scene

            clf;
            obj.simulate(simulate_filename);
            obj.save_figure([48 267 676 229],'samplex2.png')
        end

        function check_dE_from_truth (obj)
            % Absolute color errors in comparison with the ground truth
            % The left chart shows the per-patch color difference between the subject device and the ground truth. The right chart shows the bloxplot. Statistics (mean, std, min, median, and max) are provided in the titles.
            clf
            h_dE = gcf;
            ha1 = subplot(1,2,1);
            ha2 = subplot(1,2,2);
            obj.evaluate_dE_from_truth(ha1,ha2);
            set(h_dE,'Units','normalized','Position',[0 0 0.5 0.25]);
            obj.save_figure([48 267 676 229],'bar_de.png')
        end

        function check_order_axes (obj,im,ha)
            axes(ha);
            imshow(im);
        end

        function check_order (obj)
            % Order in lightness, chroma, and hue - 1D view

            clf
            hg = gcf;
            set(hg,'Visible','on')

            td = tiledlayout(4,1,'TileSpacing','tight');
            ha1 = nexttile(td,1); obj.check_order_axes(obj.im12_l6,ha1); title('Lightness order of gray patches')
            ha2 = nexttile(td,2); obj.check_order_axes(obj.im12_l,ha2); title('Lightness order of chromatic patches')
            ha3 = nexttile(td,3); obj.check_order_axes(obj.im12_h,ha3); title('Hue order of chromatic patches')
            ha4 = nexttile(td,4); obj.check_order_axes(obj.im12_c,ha4); title('Chroma order of chromatic patches')

            obj.save_figure([173 120 1731 1153],'lch_1d.png')
        end

        function check_linearity (obj)
            % Order in lightness, chroma, and hue - 2D view

            clf
            h_order = gcf();

            ha1 = subplot(2,2,1);
            obj.evaluate_order_lightness_6(ha1); title('Lightness order of gray patches')

            ha2 = subplot(2,2,2);
            obj.evaluate_order_lightness(ha2); title('Lightness order of chromatic patches')

            ha3 = subplot(2,2,3);
            obj.evaluate_order_hue(ha3); title('Hue order of chromatic patches')

            ha4 = subplot(2,2,4);
            obj.evaluate_order_chroma(ha4); title('Chroma order of chromatic patches')

            obj.save_figure([436 321 834 777],'lch_2d.png')
        end

        function check_color_transfer (obj)
            % Three-dimensional color transfer

            clf

            hg = gcf;
            set(hg,'Visible','on')

            td = tiledlayout(3,2,'TileSpacing','tight');

            ha1 = nexttile(td,1); obj.show3dquiver_select(obj.chromatic_patches,ha1); view([12 32]); axis equal; title('Chromatic patches')
            ha2 = nexttile(td,2); obj.show2dquiver_polar_select(obj.chromatic_patches,ha2); view([0 90]); axis([-100 100 -100 100 0 100]); title('Transfer in hue and chroma')
            ha3 = nexttile(td,3); obj.show3dquiver_select(obj.chromatic_patches,ha3); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('Transfer in lightness and a*')
            ha4 = nexttile(td,4); obj.show3dquiver_select(obj.chromatic_patches,ha4); view([90 0]); axis square; axis([-100 100 -100 100 0 100]); title('Transfer in lightness and b*')
            ha5 = nexttile(td,5); obj.show3dquiver_select(obj.grayscale_patches,ha5); view([0 0]); axis square; axis([-100 100 -100 100 0 100]); title('Lightness of gray patches')
            ha6 = nexttile(td,6); obj.show2dquiver_polar_select(obj.grayscale_patches,ha6); view([0 90]); axis square; axis([-30 30 -30 30 0 100]); title('Chromaticity of gray patches')

            % capture the figure
            obj.save_figure([327 131 958 1075],'labx6.png')
        end

        function check_cce (obj)
            % Preservation of Color Contrast between Patches

            clf
            h_cce = gcf;
            ha1 = subplot(1,3,1); obj.evaluate_CCE(1,ha1);
            ha2 = subplot(1,3,2); obj.evaluate_CCE(2,ha2);
            ha3 = subplot(1,3,3); obj.evaluate_CCE(3,ha3);
            set(h_cce,'Units','normalized','Position',[0 0 1 1]);

            obj.save_figure([11 529 1227 556],'ccex3.png')
        end

    end

end


