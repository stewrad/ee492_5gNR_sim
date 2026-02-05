classdef HARQEntity < handle
%   HQE = HARQEntity(HIDSEQ,RVSEQ,NTB) creates a HARQ entity object, HQE,
%   to manage a set of parallel HARQ processes for a single UE. 
%   HIDSEQ defines the fixed sequence of HARQ process IDs scheduling the
%   stop-and-wait protocol. Optional RVSEQ specifies the RV sequence used
%   for the initial transport block transmissions and any subsequent
%   retransmissions (default is RV=0, with no retransmissions). Optional
%   NTB specifies the number of transport blocks to manage for each process
%   (default is 1). A HARQEntity object is a handle object representing
%   one HARQ entity per DL-SCH/UL-SCH transport channel per UE MAC entity.
%
%   HARQ information for DL-SCH or for UL-SCH transmissions consists of,
%   new data indicator (NDI), transport block size (TBS), redundancy 
%   version (RV), and HARQ process ID. A HARQEntity object stores the HARQ
%   information for each of a set of parallel HARQ processes. Only one
%   process is active for data transmission at a time. The object steps
%   through these processes using a fixed sequence and the object 
%   properties contain the HARQ information for the currently active 
%   process. The TBS and CRC error of the associated data transmission 
%   update the current HARQ process state.
% 
%   HARQEntity properties (Read-only):
% 
%   HARQProcessID        - Current HARQ process ID
%   TransportBlockSize   - Current TBS for process, per codeword
%   TransmissionNumber   - Current transmission number in RV sequence, per codeword (0-based)
%   RedundancyVersion    - Current redundancy version, per codeword
%   CBGTI                - Current code block group transmission information, per codeword
%   NewData              - Is this new transport data for the process i.e. is this the start of a new RV sequence?
%   SequenceTimeout      - Is this a new data because of timeout i.e. last sequence ended without a successful transmission?
%
%   HARQEntity methods:
%
%   HARQEntity           - Create a HARQEntity object
%   updateProcess        - Update current HARQ process with data transmission information (TBS, CRC error, bit capacity)
%   advanceToNextProcess - Advance entity to next HARQ process in the sequence
%   updateAndAdvance     - Update current HARQ process and advance to the next    

%   Copyright 2021-2025 The MathWorks, Inc.
   
    properties (SetAccess = private)
        HARQProcessID;      % Current HARQ process ID
        TransportBlockSize  % Current TBS, per codeword
        TransmissionNumber; % Current transmission number in RV sequence (0-based)
        RedundancyVersion;  % Current redundancy version
        CBGTI;              % Current code block group transmission information, per codeword

        NewData;            % Is this the start of a new RV sequence?
        SequenceTimeout;    % Is this a new RV sequence because of timeout i.e. did the last sequence end without a successful transmission?
    end
       

    properties (Access = private)
        HarqProcessIDOrder;        % HARQ process ID sequence order to be used
        HarqRVSequence;            % RV sequence shared by all processes

        HarqProcessStateOrder;     % HARQ process sequence state order (maps the ID order into state array)
        HarqProcessStates;         % Array of state information for the individual processes
        HarqProcessIndex;          % Current index into the HARQ process sequence order
    end

    properties (Access = private)
        TotalBlocks;               % Total number of transport blocks sent/received across all processes
        SuccessfulBlocks;          % Total number of transport blocks received successfully across all processes
       
        TotalBits;                 % Total number of information bits sent/received across all processes
        SuccessfulBits;            % Total number of information received successfully across all processes
    end
    
    methods

        function obj = HARQEntity(processorder,rvsequence,ncw,maxncbg)
        %HARQEntity Create a HARQEntity object for a fixed process sequence
        %   HQE = HARQEntity(HIDSEQ,RVSEQ,NTB,MAXNCBG) creates a HARQ entity object, HQE,
        %   to manage a set of parallel HARQ processes for a single UE.
        %   HIDSEQ defines the fixed sequence of HARQ process IDs scheduling the
        %   stop-and-wait protocol. Optional RVSEQ specifies the RV sequence used
        %   for the initial transport block transmissions and any subsequent
        %   retransmissions (default is RV=0, with no retransmissions). Optional
        %   NTB specifies the number of transport blocks to manage for each process
        %   (default is 1). Optional MAXNCBG specifies the maximum number of code
        %   block groups to be used in the case of CBG-based transmission 
        %   (default is TB-based transmission, equivalent to MAXNCBG = 1 and a 
        %   single group of all code blocks). A HARQEntity object is a handle object
        %   representing one HARQ entity per DL-SCH/UL-SCH transport channel 
        %   per UE MAC entity.
            
            if nargin < 4
                maxncbg = 1;  % Default to TB based retransmission (a single group containing all code blocks)  
                if nargin < 3
                    ncw = 1;    % Default to managing a single codeword
                    if nargin < 2
                        rvsequence = 0;  % Default to no HARQ retransmissions (initial transmission only)
                    end
                end
            end

            % Store the common RV sequence for all processes
            obj.HarqRVSequence = rvsequence;     % Share the same RV sequence per CW and across all processes
    
            % Store the HARQ process fixed sequence order 
            obj.HarqProcessIDOrder = processorder;
            % Get the number of unique processes in the sequence 
            % and a 1-based index for each process ID, which will be 
            % used to link each ID to its individual process state
            [processids,~,obj.HarqProcessStateOrder] = unique(obj.HarqProcessIDOrder);        
    
            % Create array of HARQ process states for the set of the processes
            processstate.RVIdx = zeros(1,ncw);      % Indices into RV sequence for codewords
            processstate.Timeout = zeros(1,ncw);    % Indication of previous RV sequence timeout
            processstate.TBS = zeros(1,ncw);        % Transport block sizes
            processstate.CBGTI = ones(maxncbg,ncw); % Initialise the CBGTI 'bitmap' to all CBG sent (initial transmission)
            obj.HarqProcessStates = repmat(processstate,numel(processids),1);  % Create array
    
            % Initialise current process index and load its state into the public properties
            obj.HarqProcessIndex = 1;
            dereferenceProcessIndex(obj);
    
            % Clear block transmission counters
            obj.TotalBlocks = zeros(1,ncw);
            obj.SuccessfulBlocks = zeros(1,ncw);
            obj.TotalBits = zeros(1,ncw);
            obj.SuccessfulBits = zeros(1,ncw);

        end   
        
        function rep = updateProcess(obj,txerror,tbs,g,varargin)
        %updateProcess Update current HARQ process on data transmission
        %   TR = updateProcess(HE,BLKERR,TBS,G,CBGERR) updates the current HARQ process 
        %   with the per-transport block CRC error, BLKERR, the transport block size, TBS, 
        %   and the bit capacity, G, of the associated data transmission. The optional 
        %   CBGERR 1 or 2 column input is the per-transport block CBG errors for the case
        %   of CBG-based transmission (when the object was constructed with MAXNCBG > 1). 
        %   The TR output is a text update report of the result of the data transmission
        %   on the process state.
        %  
        %   Example:
        %   harqidseq = 0:15;
        %   rvseq = [0 2 3 1];
        %   harqent1 = HARQEntity(harqidseq,rvseq)
        %   urep1 = updateProcess(harqent1,1,100,300);
        %   urep1
        %   harqent1
        % 
        %   harqent2 = HARQEntity(harqidseq,rvseq,2)
        %   urep2 = updateProcess(harqent2,[1 0],[100 200],[300 600]);
        %   urep2
        %   harqent2
        %
        %   See also updateAndAdvance.
    
            % Local function for adjusting the input parameter sizes to be consistent with NCW
            ncw = numel(obj.NewData);
            function out = resizeondemand(in,mdata)   

                if nargin == 1 || mdata == 0
                   % Treat data as a row vector
                   % If a matrix data is not supported in this call
                   % then force data into a scalar/row
                   in = in(:).';
                end
             
                if iscolumn(in)  % Scalar/column
                    out = repmat(in,1,ncw); % Single column in or scalar in, but expand to 'ncw' (1 or 2) in length (output is a vector)
                else
                    % Working with a row or multi-column matrix now
                    % so limit the number of columns to match ncw
                    out = in(:,1:ncw);      % Vector in, but limit to 'ncw' (1 or 2) in length (output is the same as the input)
                end
            end
                
            % Scalar expand or contract the input parameter dimensions for NCW
            txerror = resizeondemand(txerror);
            tbs = resizeondemand(tbs);
            g = resizeondemand(g);
            
            % Prepare information related the CBG in use, if in CBG-based transmission mode
            ncb = zeros(1,ncw);
            maxncbg = size(obj.CBGTI,1);
            if ~isempty(varargin)
                cbgerr = resizeondemand(varargin{1},1);     
                if size(cbgerr,1) ~= maxncbg
                    error('CBGERR input must have the same number of rows (%d) as the CBGTI property (%d).',size(cbgerr,1),maxncbg);
                end
                % Get the number of code blocks for use later
                for cw=1:ncw
                    trinfo = nrDLSCHInfo(tbs(cw),1023/1024);
                    ncb(cw) = trinfo.C;
                end
            else
                % If no CBG error info was provided then apply the TB error 
                cbgerr = txerror & obj.CBGTI;
            end

            % Process the result of the _current_ transmission, given the
            % combination of the shared channel configuration (subset) and
            % the resulting CRC error
            if nargout
                rep = createUpdateReport(obj,txerror,tbs,g,cbgerr,ncb);
            end

            % Create a text summary of what happened for the transmission event
            % for the current process 

            % Update current HARQ process information (this updates the RV
            % depending on CRC pass or fail in the previous transmission for
            % this HARQ process)

            % Get index for current process in the state array 
            stateidx = obj.HarqProcessStateOrder(obj.HarqProcessIndex);
     
            %
            % Update the primary HARQ TB state, given the error state, ready for the _next_ transmission.
            %

            % Capture the TBS values of any new data transmissions that occurred
            obj.HarqProcessStates(stateidx).TBS(obj.NewData) = tbs(obj.NewData);

            % Check that the TBS values of any retransmissions are the same as the initial transmission
            d = (tbs ~= obj.HarqProcessStates(stateidx).TBS) & ~obj.NewData;
            if any(d)
                warning('For HARQ process %d, transport block sizes of a retransmission (%s) changed from the initial transmission (%s).',...
                     obj.HARQProcessID,...
                     join(string(tbs(d)),','),...
                     join(string(obj.HarqProcessStates(stateidx).TBS(d)),','));
            end

            % Get RV indices for current process in the RV sequence
            rvidx = obj.HarqProcessStates(stateidx).RVIdx;

            % Update process information 
            % If TB error, then advance and RV index and check for timeout, otherwise no TB error so set index back to 0
            inerror = txerror ~= zeros(size(rvidx));            % Allow 'scalar expansion' of error across codewords
            rvidx(~inerror) = 0;                                % Reset sequence indices if no error
            rvidx = rvidx + inerror;                            % Or increment if error
            timeout = rvidx == length(obj.HarqRVSequence);      % Test for a sequence timeout when we increment
            rvidx(timeout) = 0;                                 % And reset any effected by timeout

            % 
            % Update the CBGTI state, given the error state, ready for the _next_ transmission.
            % 

            % Only process the CBG information if using more than one group, otherwise in standard TB-based transmission mode
            if maxncbg > 1
                cbgti = obj.CBGTI;
                % Check contents of cbgerr for consistency for the current CBGTI
                pi = xor(cbgti,cbgerr);
                tii = find(pi);
                noteebutee = tii(cbgti(tii)==0);  % The CBGTI send bit was off, but the associated CBG error bit was on 
                if ~isempty(noteebutee)
                    [r,c]=ind2sub(size(cbgti),noteebutee);
                    etxt = sprintf("(cw,cbg)=%s",join(compose("(%d,%d)",[c,r]),','));
                    warning('Errors reported in CBGERR that are not enabled by the current CBGTI and TBS of the HARQ process: %s. These will be ignored in the HARQ process updating.',etxt); 
                end

                % Retransmit all the CBG that were marked for send AND then received in error
                nextcbgti = cbgti & cbgerr;
                obj.HarqProcessStates(stateidx).CBGTI = nextcbgti;

                % However, the new CBGTI state per TB/CW, based on the CBG error feedback alone, may require resetting under the following conditions,
                % - If the RV sequence has restarted (i.e. no TB error) then all CBG should immediately be marked for send         
                % - If the RV sequence has not restarted (i.e. a TB error), but the updated CBGTI now has no CBG marked for send (i.e. all zeros), then all CBG should be marked for send

                % Detect whether ALL the updated, TBS-active CBG bits are now marked off for resend
                nocbgfailed = zeros(1,ncw);
                for cw=1:ncw
                    m = min(ncb(cw),maxncbg);           % Number of CBG in active play for this TBS
                    nocbgfailed(cw) = all(obj.HarqProcessStates(stateidx).CBGTI(1:m,cw)==0,1);
                end
                % Detect whether all active CBG resend bits are off but there was still a TB error
                tbfailbutnocbgfail = logical(txerror) & nocbgfailed;
    
                % If a TB fail but no fails at all on the CBG set then mark all CBG for full transmission again.
                % If any rvidx that are now sitting at 0 then this will be an initial transmissions, by definition,
                % so mark all CBG for transmission
                obj.HarqProcessStates(stateidx).CBGTI(:,tbfailbutnocbgfail | (rvidx==0)) = 1;
            end

            % Capture results into the process state
            obj.HarqProcessStates(stateidx).RVIdx = rvidx;
            obj.HarqProcessStates(stateidx).Timeout = timeout;

            % Assume that each method call is for a separate transmission event
            obj.SuccessfulBlocks = obj.SuccessfulBlocks + ~inerror;
            obj.TotalBlocks = obj.TotalBlocks + 1;

            obj.SuccessfulBits = obj.SuccessfulBits + tbs.*(~inerror);   % tbs
            obj.TotalBits = obj.TotalBits + tbs;                         % tbs

            % Reflect updated state back into the public properties
            dereferenceProcessIndex(obj);

        end
        
        function advanceToNextProcess(obj)
        %advanceToNextProcess Advance to next HARQ process in process ID sequence
        %   advanceToNextProcess(HE) moves the HARQ entity on to the next process in
        %   the HARQ ID sequence. The object properties are updated with the new 
        %   process HARQ information state.
        %
        %   See also updateAndAdvance.

            % Update process index in the ID sequence order
            obj.HarqProcessIndex = mod(obj.HarqProcessIndex,length(obj.HarqProcessIDOrder))+1;  % 1-based
                        
            % Reflect updated state back into the public properties
            dereferenceProcessIndex(obj);
        end

        function varargout = updateAndAdvance(obj,error,tbs,g,varargin)
        %updateAndAdvance Update current HARQ process and advance to next process
        %   TR = updateAndAdvance(HE,BLKERR,TBS,G) updates the current HARQ process 
        %   with associated data transmission, and advance to the next process
        %   in the HARQ process ID sequence. The current process is updated with
        %   the per-transport block CRC error, BLKERR, transport block size, TBS, and
        %   bit capacity, G, of the associated data transmission. TR is a text
        %   update report of the result of the data transmission on the process state.
        %
        %   See also updateProcess, advanceToNextProcess.

            [varargout{1:nargout}] = updateProcess(obj,error,tbs,g,varargin{:});
            advanceToNextProcess(obj);
        end

    end

    methods (Access = private)
        
        % Reflect current state back in the public properties
        function dereferenceProcessIndex(obj)
            
            % Load the current indexed process state into the public properties
            obj.HARQProcessID = obj.HarqProcessIDOrder(obj.HarqProcessIndex);       % Current HARQ process ID

            stateidx = obj.HarqProcessStateOrder(obj.HarqProcessIndex);
 
            obj.TransportBlockSize = obj.HarqProcessStates(stateidx).TBS;           % Current TBS
            obj.TransmissionNumber = obj.HarqProcessStates(stateidx).RVIdx;         % Current transmission number in sequence
            obj.RedundancyVersion = obj.HarqRVSequence(obj.TransmissionNumber + 1); % Current redundancy version
            obj.CBGTI = obj.HarqProcessStates(stateidx).CBGTI;                      % Current CBG transmission information

            obj.NewData = obj.TransmissionNumber == 0;                              % Is this the start of a new sequence?
            obj.SequenceTimeout = obj.HarqProcessStates(stateidx).Timeout;          % Did the last sequence end without a successful transmission
        end
        
        % Create a text report of the effect of the data transmission on the current HARQ process
        function sr = createUpdateReport(harqEntity,blkerr,tbs,g,varargin)
            
            % Display transport block CRC error information per codeword managed by current HARQ process
            icr = tbs./g;    % Instantaneous code rate for this individual transmission
            blkerr = logical(blkerr);
            
            % Leading intro part
            strparts = sprintf("HARQ Proc %d:",harqEntity.HARQProcessID);
    
            estrings = ["passed","failed"];
            for cw=1:length(harqEntity.NewData)
                % Create a report on the RV state given position in RV sequence and decoding error
    
                % Transmission number part
                if harqEntity.NewData(cw)
                    ts1 = sprintf("Initial transmission");
                else
                    ts1 = sprintf("Retransmission #%d",harqEntity.TransmissionNumber(cw));
                end
    
                % CBG related part, if more than one CBG was present in the CBGTI state i.e. not TB-based
                ts2 = '';
                cbgsent = harqEntity.CBGTI;   % Scheduled CBGTI 
                maxncbg = size(cbgsent,1);    % Assume this dimension represents the max number of CBG
                if ~isempty(varargin) && maxncbg > 1

                    % Get the CBG related optional inputs
                    cbginerror = varargin{1};
                    ncb = varargin{2};

                    % Derive the number of CBG sent/received in error
                    m = min(ncb(cw),maxncbg);                        % Number of CBG in active play for this TBS
                    cbgs = sum(cbgsent(1:m,cw));                     % Total number of CBG sent
                    cbgse = sum(cbgsent(1:m,cw)&cbginerror(1:m,cw)); % Total number of CBG sent that were received in error

                    % CBG related dimensionality, derived from the TBS/NCB
                    m1 = mod(ncb(cw),m);    % Number of CBG of size k1 (the first m1 CBG)
                    m2 = m - m1;            % number of CBG of size k2 (the next m2 CBG)
                    k1 = ceil(ncb(cw)/m);   % k1 CB in the first set of CBG
                    k2 = floor(ncb(cw)/m);  % k2 CB in the second set of CBG
                   
                    % Total number of CB associated with the set of CBG sent
                    ncbs = [k1*ones(1,m1), k2*ones(1,m2)] * cbgsent(1:m,cw);  % Number of CB sent  

                    % ICR above is full TBS vs the actual bit capacity 
                    icr(cw) = icr(cw)*(ncbs/ncb(cw));  % Scale the full TB code rate by the ratio of CB sent to CB in full TB
         
                    % Supply additional information in the event of a discrepancy between
                    % the TB error and the CBG errors 
                    extxt = [];
                    if logical(cbgse) ~= blkerr(cw)
                        if blkerr(cw)
                            % A TB error but 0 CBG errors reported
                            extxt = "however TB in error ";
                        else
                            % No TB error but some CBG errors reported
                            extxt = "however TB passed ";
                        end
                    end

                    ts2 = sprintf("(%d of %d CBG sent failed %s(NCBG=%d,NCB=%d))",cbgse,cbgs,extxt,m,ncb(cw));
                end

                % Transmission parameters part
                ts3 = sprintf("(TBS=%d,RV=%d,CR=%f)",tbs(cw),harqEntity.RedundancyVersion(cw),icr(cw));   % For existing info, would need csn and icr

                % Add codeword report to list of string parts
                strparts(end+1) = sprintf("CW%d: %s %s %s %s.",cw-1,ts1,estrings{1+blkerr(cw)},ts2,ts3); %#ok<AGROW> 
            end
    
            % Combine all string parts
            sr = join(strparts,' ');
        
        end   
    end
   
end