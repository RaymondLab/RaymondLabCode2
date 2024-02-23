;-----------------------------------------------------------------------------
; Initial Parameters
;-----------------------------------------------------------------------------
            ; Set rate to 1 ms/step, 1 scaling, and 0 offset            
            SET    1.000,1,0

            ; Flags (1=running, 0=stopped)
            VAR    V1,BlockFlg=0   ;General block flag
            VAR    V2,Counter1     ;Main counter
            VAR    V3,Counter2     ;Sub-counter 1
            VAR    V4,Counter3     ;Sub-counter 2
            VAR    V5,Counter4     ;Sub-counter 3
            VAR    V6,Gap1Dur=478  ;(Gap1Dur*rate/5)-2
            VAR    V7,FlashDur=39  ;(FlashDur*rate/5)-1
            VAR    V8,Gap3Dur=479  ;(Gap3Dur*rate/5)-1

            ; Drum parameters
            VAR    V9,DrumOff=0
            VAR    V10,DrumAmp=0
            VAR    V11,DrumNeg=0
            VAR    V12,DRampSlp=0
            VAR    V13,DStepDur=50
            VAR    V14,DPawsDur=25
            VAR    V15,DDelDur=0

            ; Chair parameters
            VAR    V16,ChairOff=0
            VAR    V17,ChairAmp=0
            VAR    V18,ChairNeg=0
            VAR    V19,CRampSlp=0
            VAR    V20,CStepDur=50
            VAR    V21,CPawsDur=25
            VAR    V22,CDelDur=0

            ; Block parameters
            VAR    V23,NumTest=3
            VAR    V24,NumTrain=5

            ; TTL parameters
            VAR    V25,OptoStrt=0
            VAR    V26,OptoDur=0
            VAR    V27,OptoInt=0
            VAR    V28,OptoCycl=0
            VAR    V29,LiteType=1
            VAR    V30,SkipPos=1
            VAR    V31,SkipNeg=1


;-----------------------------------------------------------------------------
; IDLELOOP: Sequencer loop when all signals are OFF
;-----------------------------------------------------------------------------
IDLELOOP: 'I DAC   0,DrumOff       ;Idling loop        >IDLING
            DAC    1,ChairOff      ;Repeatedly sets offset of drum and chair >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TTL1OFF: Turns off TTL channel 1 (digital output bit 8)
;-----------------------------------------------------------------------------
TTL1OFF: 'l DIGOUT [.......0]      ;Turn TTL 1 off     >=


;-----------------------------------------------------------------------------
; TTL1ON: Turns on TTL channel 1 (digital output bit 8)
;-----------------------------------------------------------------------------
TTL1ON: 'L  DIGOUT [.......1]      ;Turn TTL 1 on      >=


;-----------------------------------------------------------------------------
; RESET: Resets sequencer to initial state
;-----------------------------------------------------------------------------
RESET:  'R  DIGOUT [.......0]      ;Reset to initial state >=
            RAMP   1,ChairOff,CRampSlp ;Ramp chair voltage back to offset >=
            RAMP   0,DrumOff,DRampSlp ;Ramp drum voltage back to offset >=
RESET1:     WAITC  1,RESET1        ;Wait for chair ramp to finish >=
RESET2:     WAITC  0,RESET2        ;Wait for drum ramp to finish >=
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; DSTEP: Single positive / negative drum velocity step
;-----------------------------------------------------------------------------
DSTEP:  'D  RAMP   0,DrumAmp,DRampSlp ;Start positive ramp >=
DRAMP1:     WAITC  0,DRAMP1        ;Wait for ramp to finish >=
            MOV    Counter3,DStepDur ;Single drum step >=
DWAIT1:     DAC    0,DrumAmp       ;Set drum amplitude >=
            DAC    1,ChairOff      ;Set chair offset   >=
            DBNZ   Counter3,DWAIT1 ;Repeat until counter hits zero >=

            RAMP   0,DrumOff,DRampSlp ;Start ramp to offset >=
DRAMP2:     WAITC  0,DRAMP2        ;Wait for ramp to finish >=
            BEQ    DPawsDur,0,DSKIP1 ;If no pause, go to DSKIP1 >=
            MOV    Counter3,DPawsDur ;Set drum step pause duration >=
DWAIT2:     DAC    0,DrumOff       ;Set drum offset    >=
            DAC    1,ChairOff      ;Set chair offset   >=
            DBNZ   Counter3,DWAIT2 ;Repeat until counter hits zero >=

DSKIP1:     RAMP   0,DrumNeg,DRampSlp ;Start negative ramp >=
DRAMP3:     WAITC  0,DRAMP3        ;Wait for ramp to finish >=
            MOV    Counter3,DStepDur ;Set negative drum step duration >=
DWAIT3:     DAC    0,DrumNeg       ;Set drum amplitude >=
            DAC    1,ChairOff      ;Set chair amplitude >=
            DBNZ   Counter3,DWAIT3 ;Repeat until counter hits zero >=

            RAMP   0,DrumOff,DRampSlp ;Start ramp to offset >=
DRAMP4:     WAITC  0,DRAMP4        ;Wait for ramp to finish >=
            BEQ    DPawsDur,0,DSKIP2 ;If no pause, go to DSKIP2 >=
            MOV    Counter3,DPawsDur ;Set drum step pause duration >=
DWAIT4:     DAC    0,DrumOff       ;Set drum offset    >=
            DAC    1,ChairOff      ;Set chair offset   >=
            DBNZ   Counter3,DWAIT4 ;Repeat until counter hits zero >=
DSKIP2:     NOP    
            RETURN 


;-----------------------------------------------------------------------------
; CSTEP: Single positive / negative chair velocity step
;-----------------------------------------------------------------------------
CSTEP:  'C  RAMP   1,ChairAmp,CRampSlp ;Start positive ramp >=
CRAMP1:     WAITC  1,CRAMP1        ;Wait for ramp to finish >=
            MOV    Counter3,CStepDur ;Single chair step >=
CWAIT1:     DAC    1,ChairAmp      ;Set chair amplitude >=
            DAC    0,DrumOff       ;Set drum offset    >=
            DBNZ   Counter3,CWAIT1 ;Repeat until counter hits zero >=

            RAMP   1,ChairOff,CRampSlp ;Start ramp to offset >=
CRAMP2:     WAITC  1,CRAMP2        ;Wait for ramp to finish >=
            BEQ    CPawsDur,0,CSKIP1 ;If no pause, go to CSKIP1 >=
            MOV    Counter3,CPawsDur ;Set chair step pause duration >=
CWAIT2:     DAC    1,ChairOff      ;Set chair offset   >=
            DAC    0,DrumOff       ;Set drum offset    >=
            DBNZ   Counter3,CWAIT2 ;Repeat until counter hits zero >=

CSKIP1:     RAMP   1,ChairNeg,CRampSlp ;Start negative ramp >=
CRAMP3:     WAITC  1,CRAMP3        ;Wait for ramp to finish >=
            MOV    Counter3,CStepDur ;Set negative chair step duration >=
CWAIT3:     DAC    1,ChairNeg      ;Set chair amplitude >=
            DAC    0,DrumOff       ;Set drum amplitude >=
            DBNZ   Counter3,CWAIT3 ;Repeat until counter hits zero >=

            RAMP   1,ChairOff,CRampSlp ;Start ramp to offset >=
CRAMP4:     WAITC  1,CRAMP4        ;Wait for ramp to finish >=
            BEQ    CPawsDur,0,CSKIP2 ;If no pause, go to CSKIP2 >=
            MOV    Counter3,CPawsDur ;Set chair step pause duration >=
CWAIT4:     DAC    1,ChairOff      ;Set chair offset   >=
            DAC    0,DrumOff       ;Set drum offset    >=
            DBNZ   Counter3,CWAIT4 ;Repeat until counter hits zero >=
CSKIP2:     NOP                    ;                   >=
            RETURN 


;-----------------------------------------------------------------------------
; TRAIN: x2/x1/x0 training block
;-----------------------------------------------------------------------------
TRAIN:  'T  MOVI   BlockFlg,1      ;Start training block >TRAINING
            MOV    Counter3,NumTrain ;Set number of step periods to run >"
            BNE    LiteType,1,TRAIN0 ;Check whether light always on or opto >"
            DIGOUT [.......1]      ;Turn TTL 1 on      >"
TRAIN0:     BEQ    SkipPos,1,TSKIP1 ;Check whether positive opto or skip >"
            DIGPS  0,P,OptoStrt    ;Set pulse start time + (OptoDur/2) >"
            DIGPS  0,D,OptoDur     ;Set pulse duration >"
            DIGPS  0,C,OptoCycl    ;Number of cycles   >"
            DIGPC  0,G             ;Start pulse train  >"
TSKIP1:     RAMP   0,DrumAmp,DRampSlp ;                >"
            RAMP   1,ChairAmp,CRampSlp ;               >"
TRAIN1:     WAITC  0,TRAIN1        ;                   >"
TRAIN2:     WAITC  1,TRAIN2        ;                   >"
            DIGPS  0,I,OptoInt     ;                   >"
            MOV    Counter1,DStepDur ;                 >"
TRAIN3:     DAC    0,DrumAmp       ;                   >"
            DAC    1,ChairAmp      ;                   >"
            DBNZ   Counter1,TRAIN3 ;Run steps until counter hits zero >"

            RAMP   0,DrumOff,DRampSlp ;                >"
            RAMP   1,ChairOff,CRampSlp ;               >"
TRAIN4:     WAITC  0,TRAIN4        ;                   >"
TRAIN5:     WAITC  1,TRAIN5        ;                   >"
            BLE    DPawsDur,0,TSKIP2 ;If no pause, go to TSKIP2 >=
            MOV    Counter1,DPawsDur ;                 >"
TRAIN6:     DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   Counter1,TRAIN6 ;Run steps until counter hits zero >"


TSKIP2:     BEQ    SkipNeg,1,TSKIP3 ;                  >"
            DIGPS  0,P,OptoStrt    ;                   >"
            DIGPS  0,D,OptoDur     ;                   >"
            DIGPS  0,C,OptoCycl    ;                   >"
            DIGPC  0,G             ;                   >"
TSKIP3:     RAMP   0,DrumNeg,DRampSlp ;                >"
            RAMP   1,ChairNeg,CRampSlp ;               >"
TRAIN7:     WAITC  0,TRAIN7        ;                   >"
TRAIN8:     WAITC  1,TRAIN8        ;                   >"
            DIGPS  0,I,OptoInt     ;                   >"
            MOV    Counter1,DStepDur ;Set number of step periods to run >"
TRAIN9:     DAC    0,DrumNeg       ;                   >"
            DAC    1,ChairNeg      ;                   >"
            DBNZ   Counter1,TRAIN9 ;Run steps until counter hits zero >"

            RAMP   0,DrumOff,DRampSlp ;                >"
            RAMP   1,ChairOff,CRampSlp ;               >"
TRAIN10:    WAITC  0,TRAIN10       ;                   >"
TRAIN11:    WAITC  1,TRAIN11       ;                   >"
            BLE    DPawsDur,0,TSKIP4 ;If no pause, go to TSKIP4 >=
            MOV    Counter1,DPawsDur ;                 >"
TRAIN12:    DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   Counter1,TRAIN12 ;Run steps until counter hits zero >"
TSKIP4:     NOP    
            DBNZ   Counter3,TRAIN0 ;                   >"
            DIGOUT [.......0]      ;Ensure TTL 1 is off >=
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TEST: pre/post test VORD (chair only) block
;-----------------------------------------------------------------------------
TEST:   'P  MOVI   BlockFlg,1      ;Start Test block   >TESTING
            DIGOUT [.......0]      ;Ensure light is off >"
            MOV    Counter2,NumTest ;Set number of step periods to run >"
TEST1:      CALL   CSTEP           ;Run a chair step   >"
            DBNZ   Counter2,TEST1  ;Run steps until counter hits zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; GAP: No stimuli for Gap1Dur+Gap3Dur ms and FlashDur ms light pulse halfway
;-----------------------------------------------------------------------------
GAP:    'G  MOVI   BlockFlg,1      ;Start Gap block    >GAP BLOCK
            DIGOUT [.......0]      ;Reset to initial state >"
            RAMP   1,ChairOff,CRampSlp ;Ramp chair voltage back to offset >"
            RAMP   0,DrumOff,DRampSlp ;Ramp drum voltage back to offset >"
GAPINIT1:   WAITC  1,GAPINIT1      ;Wait for chair ramp to finish >"
GAPINIT2:   WAITC  0,GAPINIT2      ;Wait for drum ramp to finish >"
            MOV    Counter1,Gap1Dur ;Set duration of first half >"
GAP1:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   Counter1,GAP1   ;Repeat until counter hits zero >"
            BEQ    FlashDur,0,GSKIP ;Skip to GAP2 if FlashDur is zero >"
            DIGOUT [.......1]      ;Turn light on      >"
            MOV    Counter1,FlashDur ;Set duration of light pulse >"
GAP2:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   Counter1,GAP2   ;Repeat until counter hits zero >"
GSKIP:      DIGOUT [.......0]      ;Turn light off     >"
            MOV    Counter1,Gap3Dur ;Set duration of second half >"
GAP3:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   Counter1,GAP3   ;Repeat until counter hits zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP