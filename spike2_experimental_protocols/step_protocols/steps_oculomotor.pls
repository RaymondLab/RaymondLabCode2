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
            VAR    V5,Gap1Dur=478  ;(Gap1Dur*rate/5)-2
            VAR    V6,FlashDur=39  ;(FlashDur*rate/5)-1
            VAR    V7,Gap3Dur=479  ;(Gap3Dur*rate/5)-1

            ; Drum parameters
            VAR    V8,DrumOff=0
            VAR    V9,DrumAmp=0
            VAR    V10,DrumNeg=0
            VAR    V11,DrumFrq=0
            VAR    V12,DrumPha=0

            ; Chair parameters
            VAR    V13,ChairOff=0
            VAR    V14,ChairAmp=0
            VAR    V15,ChairNeg=0
            VAR    V16,ChairFrq=0
            VAR    V17,ChairPha=0

            VAR    V18,LightPer
            VAR    V19,StepDel=100
            VAR    V20,StepDur=100
            VAR    V21,StepOff=50
            VAR    V22,Slope=0


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
            RAMP   1,ChairOff,Slope ;Ramp chair voltage back to offset >=
            RAMP   0,DrumOff,Slope ;Ramp drum voltage back to offset >=
RESET1:     WAITC  1,RESET1        ;Wait for chair ramp to finish >=
RESET2:     WAITC  0,RESET2        ;Wait for drum ramp to finish >=
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; DSTEP: Single positive / negative drum velocity step
;-----------------------------------------------------------------------------
DSTEP:  'D  BLE    StepDel,0,DSKIP1 ;Single drum step  >=
            MOV    Counter1,StepDel ;                  >=
DWAIT1:     DAC    0,DrumOff       ;IDLELOOP during initial delay >=
            DAC    1,ChairOff      ;                   >=
            DBNZ   Counter1,DWAIT1 ;Repeat until counter hits zero >=
DSKIP1:     MOV    Counter1,StepDur ;Set positive drum step duration >=
            RAMP   0,DrumAmp,Slope ;Start positive ramp >=
DRAMP1:     WAITC  0,DRAMP1        ;Wait for ramp to finish >=
DWAIT2:     DAC    0,DrumAmp       ;Set drum amplitude >=
            DAC    1,ChairOff      ;Set chair offset   >=
            DBNZ   Counter1,DWAIT2 ;Repeat until counter hits zero >=
            BLE    StepOff,0,DSKIP2 ;Single drum step  >=
            MOV    Counter1,StepOff ;Set drum step pause duration >=
            RAMP   0,DrumOff,Slope ;Start ramp to offset >=
DRAMP2:     WAITC  0,DRAMP2        ;Wait for ramp to finish >=
DWAIT3:     DAC    0,DrumOff       ;IDLELOOP during step pause >=
            DAC    1,ChairOff      ;                   >=
            DBNZ   Counter1,DWAIT3 ;Repeat until counter hits zero >=
DSKIP2:            MOV    Counter1,StepDur ;Set negative drum step duration >=
            RAMP   0,DrumNeg,Slope ;Start negative ramp >=
DRAMP3:     WAITC  0,DRAMP3        ;Wait for ramp to finish >=
DWAIT4:     DAC    0,DrumNeg       ;Set drum amplitude >=
            DAC    1,ChairOff      ;Set chair offset   >=
            DBNZ   Counter1,DWAIT4 ;Repeat until counter hits zero >=
            RAMP   0,DrumOff,Slope ;Start ramp to offset >=
DRAMP4:     WAITC  0,DRAMP4        ;Wait for ramp to finish >=
            RETURN 


;-----------------------------------------------------------------------------
; TRAIN: Outputs 2 TTL (channel 1) light pulses with some period
;-----------------------------------------------------------------------------
TRAIN:  'T  MOVI   BlockFlg,1      ;Oculomotor conditioning block >TRAINING
            DIGOUT [.......0]      ;Ensure light is off >"
            DIGPS  0,S,LightPer    ;Set square wave with period of T >"
            DIGPS  0,C,3           ;Run for N periods  >"
            DIGPC  0,G             ;Start TTL train    >"
TTL1:       DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DIGPBR 0,C,TTL1        ;                   >"
            CALL   DSTEP           ;                   >"
TTL2:       DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DIGPBR 0,C,TTL2        ;                   >"
            CALL   DSTEP           ;                   >"
TTL3:       DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DIGPBR 0,C,TTL3        ;                   >"
            CALL   DSTEP           ;                   >"
TTL4:       DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DIGPBR 0,S,TTL4        ;                   >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TEST: Outputs 3 TTL (channel 1) light pulses with some period
;-----------------------------------------------------------------------------
TEST:   'P  MOVI   BlockFlg,1      ;Test block         >TESTING
            DIGOUT [.......0]      ;Ensure light is off >"
            DIGPS  0,S,LightPer    ;Set square wave with period of T >"
            DIGPS  0,C,3           ;Run for N periods  >"
            DIGPC  0,G             ;Start TTL train    >"
TEST1:      DAC    0,DrumOff       ;                   >"
            DAC    1,ChairOff      ;                   >"
            DIGPBR 0,S,TEST1       ;IDLELOOP until pulses are finished >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; GAP: No stimuli for Gap1Dur+Gap3Dur ms and FlashDur ms light pulse halfway
;-----------------------------------------------------------------------------
GAP:    'G  MOVI   BlockFlg,1      ;Start Gap block    >GAP BLOCK
            DIGOUT [.......0]      ;Ensure light is off >"
            RAMP   1,ChairOff,Slope ;Ramp chair voltage back to offset >"
            RAMP   0,DrumOff,Slope ;Ramp drum voltage back to offset >"
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