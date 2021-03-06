clear all; clc; close all; % clean up
addpath(genpath('./_src'));
addpath(genpath('./_matStore'));

vdt = [1.0e-5 1.0e-6 5.0e-6];

for test_dt_ind = 1:length(vdt)

    for test_idx = 1:30

        test_cases = [500 1000 2000 3000 4000 5000 10000 25000 50000 100000];

        for stability_idx  = 1:length(test_cases)


            % Sim Parameters
            dt     = vdt(test_dt_ind);        % s


            if stability_idx <= 6
                nSpins = test_cases(stability_idx);
                nSet   = 1;
            else
                nSpins = round(test_cases(stability_idx)/32);
                nSet   = 32;
            end


            % Seq Parameters 
            pulsedurE = 0.001; %s
            pulsedurR = 0.001; %s
            TE     = 20.0e-3;       % s
            ldelta = 10e-3;        % s was 0.011
            sdelta = 5.0e-3;        % s was 0.002
            bVal   = [0,100,200,300,500,600,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,10000,12500,15000,20000]; %s/mm^2
            % bVal = 1000;%s/mm^2

            % Tissue Paramaters
            T1     = 1.5;           % s  (not used)
            T2     = 1.0;           % s  (not used)
            D      = 3.0e-9;        % m^2/s
            % Radii = [2.5e-6,3.5e-6,4.0e-6,5.0e-6,6.0e-6,8.0e-6];        % m
            Radii = [3.5e-6,4.0e-6,5.0e-6,8.0e-6];        % s/m

            %remove radius loop and do only one at a time and check

            %%
            res = zeros(nSet,length(bVal),length(Radii));
            for iSet = 1:nSet 
                fprintf('\n------------Set %d/%d ------------\n', iSet,nSet);
                for rIdx = 1:length(Radii)
                    fprintf('Radius %d/%d\n', rIdx,length(Radii));
                    %% %%%%%%%%%%%%% GET RND WALKS %%%%%%%%%%%%% 
                    spinCor = getRndWalks(D,Radii(rIdx),nSpins,(round(TE/dt) + round(pulsedurE/dt/2)),dt);

                    for bIdx = 1:length(bVal)
                        %% %%%%%%%%%%%%% GET SEQUENCE %%%%%%%%%%%%%%%
                        Sequence_Struct = getDifSeq(pulsedurE, pulsedurR, TE, sdelta, ldelta, bVal(bIdx)/1e-6, dt);
                        %check b-result

                        %% %%%%%%%%%%%%% CHECK SEQUENCE %%%%%%%%%%%%%
                        %dispSequence(Sequence_Struct);

                        %% %%%%%%%%%%%%% Perform sequence %%%%%%%%%%%
                        res(iSet,bIdx,rIdx) = simulateMRISequence(Sequence_Struct,T1,T2,spinCor,dt);

                    end
                end
            end

            %% %%%%%%%%%%%%% Final Plotting %%%%%%%%% 
            final_res = squeeze(sum(res,1)./nSet);
%             plot_signal_vs_b_results(final_res, bVal);

%             save(sprintf("./_matStoreBounce/averageRes_%d_%d_dt%d",stability_idx,test_idx,test_dt_ind),'final_res');
        end
    end
end
