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
            VAR    V11,DrumFrq=0
            VAR    V12,DrumPhs=0

            ; Chair parameters
            VAR    V13,ChairOff=0
            VAR    V14,ChairAmp=0
            VAR    V15,ChairFrq=0
            VAR    V16,ChairPhs=0
            VAR    V17,ChairAng=0

            VAR    V18,NDelays=0
            VAR    V19,PulDur=0
            VAR    V20,PulInt=0
            VAR    V21,PulNCycl=0
            VAR    V22,CycleInt=0
            VAR    V23,TTLNum1=0
            VAR    V24,TTLNum2=0

            VAR    V25,TTLNum=0
            VAR    V26,Idx1=-1
            VAR    V27,Idx2=4999
            VAR    V28,OptoDel1=-99
            VAR    V29,OptoDel2=-99

            TABSZ  10001


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
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TTL1ON: Turns on TTL channel 1 (digital output bit 8)
;-----------------------------------------------------------------------------
TTL1ON: 'L  DIGOUT [.......1]      ;Turn TTL 1 on      >=
            JUMP   IDLELOOP

;-----------------------------------------------------------------------------
; RESET: Resets sequencer to initial state
;-----------------------------------------------------------------------------
RESET:  'R  DIGOUT [.......0]      ;Reset to initial state >=
            RATE   0,0             ;Stop sine on drum  >=
            RATE   1,0             ;Stop sine on chair >=
            JUMP   IDLELOOP        ;Return to idle loop


;-----------------------------------------------------------------------------
; DRUM SINE: Turn sinusoidal drum stimulus on/off
;-----------------------------------------------------------------------------
DSINEON: 'D SZ     0,DrumAmp       ;Start drum cosine
            OFFSET 0,DrumOff       ;Set cosine offset
            PHASE  0,DrumPhs       ;Set cosine relative phase
            ANGLE  0,0             ;Set cosine angle
            RATE   0,DrumFrq       ;Set cosine frequency
DSINE1:     WAITC  0,DSINE1
            OFFSET 0,DrumOff
            JUMP   DSINE1          ;Set cosine loop    > Sine running

DSINEOFF: 'd CLRC  0               ;Stop drum sine at phase 0
DSINE2:     OFFSET 0,DrumOff       ;Adjust cosine offset
            WAITC  0,DSINE2        ;Wait for end of cycle
            RATE   0,0,IDLELOOP    ;Stop drum cosine then idle

;-----------------------------------------------------------------------------
; CHAIR SINE: Turn sinusoidal chair stimulus on/off
;-----------------------------------------------------------------------------
CSINEON: 'C SZ     1,ChairAmp      ;Start chair cosine
            OFFSET 1,ChairOff      ;Set cosine offset
            PHASE  1,ChairPhs      ;Set cosine relative phase
            ANGLE  1,0             ;Set cosine angle
            RATE   1,ChairFrq      ;Set cosine frequency
CSINE1:     WAITC  1,CSINE1
            OFFSET 1,ChairOff
            JUMP   CSINE1          ;Set cosine loop    > Sine running

CSINEOFF: 'c CLRC  1               ;Stop chair sine at phase 0
CSINE2:     OFFSET 1,ChairOff      ;Adjust cosine offset
            WAITC  1,CSINE2        ;Wait for end of cycle
            RATE   1,0,IDLELOOP    ;Stop chair cosine

;-----------------------------------------------------------------------------
; SINE: Turn sinusoidal drum and chair stimulus on/off
;-----------------------------------------------------------------------------
SINEON: 'S  SZ     0,DrumAmp       ;Start cosine
            SZ     1,ChairAmp      ;Set cosine amplitude
            OFFSET 0,DrumOff       ;Set cosine offset
            OFFSET 1,ChairOff
            PHASE  0,DrumPhs       ;Set cosine relative phase
            PHASE  1,ChairPhs
            ANGLE  0,0             ;Set cosine angle
            ANGLE  1,0
            RATE   0,DrumFrq       ;Set cosine frequency
            RATE   1,ChairFrq
SINE1:      WAITC  1,SINE1
            OFFSET 0,DrumOff
            OFFSET 1,ChairOff
            JUMP   SINE1           ;Set cosine loop    > Sine running

SINEOFF: 's CLRC   0               ;Stop sines at phase 0
            CLRC   1
SINE2:      OFFSET 0,DrumOff       ;Adjust cosine offsets
            OFFSET 1,ChairOff
            WAITC  1,SINEOFF       ;Wait for end of chair cycle
            RATE   0,0             ;Stop both cosines
            RATE   1,0
            JUMP   IDLELOOP        ;Return to idle loop

;-----------------------------------------------------------------------------
; OPTO SINE: Turn sinusoidal chair + opto stimulus on/off
;-----------------------------------------------------------------------------
OPTOON: 'O  DIGOUT [.......0]      ;Start opto stim trials
            SZ     1,ChairAmp      ;Set cosine amplitude
            OFFSET 1,ChairOff      ;Set cosine offset
            PHASE  1,ChairPhs      ;Set cosine relative phase
            ANGLE  1,ChairAng      ;Set cosine angle
            RATE   1,ChairFrq      ;Set cosine frequency

OPTO0:      MOV    Counter1,CycleInt ;Set cycle interval between stims

OPTO1:      WAITC  1,OPTO1         ;Loop until phase == ChairPhs
            OFFSET 1,ChairOff      ;Adjust cosine offset to prevent drift
            DBNZ   Counter1,OPTO1  ;Repeat until cycle counter hits zero
            BGE    Idx1,NDelays,OPTOOFF ;Branch if index exceeds number of trials
            ADDI   Idx2,1          ;Increment Idx2 by 1
            TABLD  OptoDel2,[Idx2] ;Get OptoDel2 from table at index Idx2
            ADDI   Idx1,1          ;Increment Idx1 by 1
            TABLD  OptoDel1,[Idx1] ;Get OptoDel1 from table at index Idx1
            BLT    OptoDel1,0,OPTO2 ;Skip stim if OptoDel1 < 0 ms
            MOV    TTLNum,TTLNum1  ;Otherwise, set TTLNum
            DELAY  OptoDel1        ;Delay stim by OptoDel1 number of ms
            CALL   OPTOEQ          ;Execute opto stim

OPTO2:      BLT    OptoDel2,2,OPTO0 ;Skip stim if OptoDel2 < 2 ms
            MOV    TTLNum,TTLNum2  ;Otherwise, set TTLNum
            DELAY  OptoDel2        ;Delay by some number of ms
            CALL   OPTOEQ          ;Execute opto stim
            JUMP   OPTO0           ;Set opto-cosine loop > OPTO RUNNING

OPTOEQ:     BEQ    TTLNum,2,OPTOR
OPTOL:      DIGPS  1,P,PulInt      ;Pulse every "PulInt" ms
            DIGPS  1,D,PulDur      ;Pulse has duration of "PulDur" ms
            DIGPS  1,C,PulNCycl    ;Set number of pulses in train
            DIGPC  1,G             ;Start train
            RETURN 
OPTOR:      DIGPS  2,P,PulInt      ;Pulse every "PulInt" ms
            DIGPS  2,D,PulDur      ;Pulse has duration of "PulDur" ms
            DIGPS  2,C,PulNCycl    ;Set number of pulses in train
            DIGPC  2,G             ;Start train
            RETURN 

OPTOOFF: 'o CLRC   1               ;Stop cosine DBNZ   Counter1,OPTOOFF ;Repeat until counter hits zero
OPTO3:      WAITC  1,OPTO3         ;Wait for end of cycle
            RATE   1,0,IDLELOOP    ;Stop chair cosine then idle


;-----------------------------------------------------------------------------
; WAVE: Custom waveform stimuli (if provided)
;-----------------------------------------------------------------------------
;WAVE:   'W  WAVEGO w,TW
;            WAVEST T