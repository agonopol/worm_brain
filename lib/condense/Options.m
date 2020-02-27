classdef Options
    properties
        destination = 'results/contract'
        fastStop = false;
        verbosityLevel = 1;
        indentationLevel = 0;
        numPrincipalComponents = 50;
        numDiffusionSteps = 1;
        maxNumContractionSteps = 200;
        kKNN = 20;
        initialSigma = Inf;
        nystroemN = 200;
        plotAnimation = true;
        inertia = 0.5;
        expertLabels = [];
        expId = '';
        emitRuntimeBreakdown = false;
        maxClusters = 30;
        emitCondensationSequence = false;
        mergeEpsilonClusters = true;
        thresholdControlSigma = 0.5;
        clusterAssignmentMethod;
        controlSigmaMethod;
        frequencyMergingEpsilonClusters;
        modeCalcDistances = 'normal';
        epsilonClusterIdentificationMethod = 'dynamicSigmaFraction';
        prefixFileNames = '';
        phateEmbedding = false;
        labelfn = @(clusters, labels) labels;
        sizefn = @(clusters, labels) ones(length(labels), 1);
        plotfn = @plotstate;
    end
    methods
        function obj = Options(varargin)
            mode = 'classic';
            for i=1:length(varargin)-1
                if (strcmp(varargin{i}, 'mode'))
                    mode = varargin{i+1};
                end
            end
            
            switch (mode)
                case 'classic'
                    obj.clusterAssignmentMethod = 'none';
                    obj.controlSigmaMethod = 'nuclearNormStabilization';
                    obj.frequencyMergingEpsilonClusters = 'uponMetastability';
                    
                case 'fast'
                    obj.clusterAssignmentMethod = 'none';
                    obj.controlSigmaMethod = 'movementStabilization';
                    obj.frequencyMergingEpsilonClusters = 'always';
            end
        end
        function str = asString(obj)
            str = [obj.destination ...
                   '/' ...
                   num2str(obj.numPrincipalComponents) '_' ...
                   num2str(obj.numDiffusionSteps) '_' ...
                   num2str(obj.kKNN) '_' ...
                   num2str(obj.initialSigma) '_' ...
                   num2str(obj.nystroemN) '_' ...
                   num2str(obj.inertia) '_' ...
                   num2str(obj.mergeEpsilonClusters) '_' ...
                   obj.clusterAssignmentMethod '_' ...
                   obj.controlSigmaMethod '_' ...
                   obj.frequencyMergingEpsilonClusters '_' ...
                   obj.epsilonClusterIdentificationMethod '_' ...
                   obj.modeCalcDistances ...
                  ];
        end
    end
end
