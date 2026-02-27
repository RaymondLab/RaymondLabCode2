;-----------------------------------------------------------------------------
; Initial Parameters
;-----------------------------------------------------------------------------
            ; Set rate to 1 ms/sequencer step, 1 scaling, and 0 offset
            SET    1.000,1,0

            ; --- Core Sequencer Variables ---
            ; These variables should be the same across ALL protocols
            ; Note V1-V256 is available, but DO NOT use V56-V64, V255, V256

            ; Block Flag V1 (tells Spike2 if a block is actively running)
            VAR    V1,BlockFlg=0   ;1 = running, 0 = not running

            ; Drum/Target parameters (V2 - V9)
            VAR    V2,DrumOff=0    ;Drum offset
            VAR    V3,DrumAmp=0    ;Drum amplitude
            VAR    V4,DrumFrq=0    ;Drum frequency
            VAR    V5,DrumPhs=0    ;Drum phase
            VAR    V6,DrumAng=0    ;Drum angle
            VAR    V7,NdrumP=0     ;Number of drum TEST block cycles
            VAR    V8,NdrumT=0     ;Number of drum TRAIN block cycles
            VAR    V9,DrumCtr=0    ;Cycle counter for drum

            ; Chair/Head parameters (V10 - V17)
            VAR    V10,ChairOff=0  ;Chair offset
            VAR    V11,ChairAmp=0  ;Chair amplitude
            VAR    V12,ChairFrq=0  ;Chair frequency
            VAR    V13,ChairPhs=0  ;Chair phase
            VAR    V14,ChairAng=0  ;Chair angle
            VAR    V15,NchairP=0   ;Number of chair TEST block cycles
            VAR    V16,NchairT=0   ;Number of chair TRAIN block cycles
            VAR    V17,ChairCtr=0  ;Cycle counter for chair

            ; Gap block timing parameters (V18 - V21)
            VAR    V18,Gap1Dur=478 ;(Gap1Dur*rate/5)-2
            VAR    V19,FlashDur=39 ;(FlashDur*rate/5)-1
            VAR    V20,Gap3Dur=479 ;(Gap3Dur*rate/5)-1
            VAR    V21,GapCtr=0    ;Counter for gap block

            ; Variables for storing dynamic drum/chair amp values
            VAR    V100,DrumTmp=0  ;Variable drum amplitude
            VAR    V101,ChairTmp=0 ;Variable drum amplitude

            ; --- Protocol-specific Sequencer Variables ---

            VAR    V22,NDelays=0
            VAR    V23,PulDur=0
            VAR    V24,PulInt=0
            VAR    V25,PulNCycl=0
            VAR    V26,EveryNth=0
            VAR    V27,Unused27=0  ;(was HalfCycl, no longer needed)
            VAR    V28,DelIdx=0
            VAR    V29,DelVal1=-99
            VAR    V30,DelVal2=-99
            VAR    V31,DelCh1=-99
            VAR    V32,StimNum=0
            VAR    V33,CyclCtr=0

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
            MOVI   DrumTmp,0       ;Reset drum amplitude >=
            MOVI   ChairTmp,0      ;Reset chair amplitude >=
            RATE   0,0             ;Stop sine on drum  >=
            RATE   1,0             ;Stop sine on chair >=
            JUMP   IDLELOOP        ;Return to idle loop


;-----------------------------------------------------------------------------
; DRUM SINE: Turn sinusoidal drum only stimulus on/off
;-----------------------------------------------------------------------------
DSINEON: 'D MOV    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            SZ     0,DrumAmp       ;Start drum cosine  >"
            OFFSET 0,DrumOff       ;Set cosine offset  >"
            PHASE  0,-90           ;Set cosine relative phase >"
            ANGLE  0,0             ;Set cosine angle   >"
            RATE   0,DrumFrq       ;Set cosine frequency
DSINE1:     WAITC  0,DSINE1        ;                   >"
            OFFSET 0,DrumOff       ;                   >"
            JUMP   DSINE1          ;Set cosine loop    > Drum running

DSINEOFF: 'd CLRC  0               ;Stop drum sine at phase 0 >"
DSINE2:     OFFSET 0,DrumOff       ;Adjust cosine offset >"
            WAITC  0,DSINE2        ;Wait for end of cycle >"
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            RATE   0,0,IDLELOOP    ;Stop drum cosine then idle


;-----------------------------------------------------------------------------
; CHAIR SINE: Turn sinusoidal chair only stimulus on/off
;-----------------------------------------------------------------------------
CSINEON: 'C MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 1,ChairOff      ;Set cosine offset  >"
            PHASE  1,-90           ;Set cosine relative phase >"
            ANGLE  1,0             ;Set cosine angle   >"
            RATE   1,ChairFrq      ;Set cosine frequency >"
CSINE1:     WAITC  1,CSINE1        ;                   >"
            OFFSET 1,ChairOff      ;                   >"
            JUMP   CSINE1          ;Set cosine loop    > Chair running

CSINEOFF: 'c CLRC  1               ;Stop chair sine at phase 0 >"
CSINE2:     OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,CSINE2        ;Wait for end of cycle >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            RATE   1,0,IDLELOOP    ;Stop chair cosine


;-----------------------------------------------------------------------------
; SINE: Turn sinusoidal drum and chair stimulus on/off
;-----------------------------------------------------------------------------
SINEON: 'S  MOV    DrumTmp,DrumAmp ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            SZ     0,DrumAmp       ;Start cosine       >"
            SZ     1,ChairAmp      ;Set cosine amplitude >"
            OFFSET 0,DrumOff       ;Set cosine offset  >"
            OFFSET 1,ChairOff
            PHASE  0,-90           ;Set cosine relative phase >"
            PHASE  1,-90
            ANGLE  0,0             ;Set cosine angle   >"
            ANGLE  1,0
            RATE   0,DrumFrq       ;Set cosine frequency >"
            RATE   1,ChairFrq
SINE1:      WAITC  1,SINE1
            OFFSET 0,DrumOff
            OFFSET 1,ChairOff
            JUMP   SINE1           ;Set cosine loop    > Sine running

SINEOFF: 's CLRC   0               ;Stop sines at phase 0 of drum >"
SINE2:      OFFSET 0,DrumOff       ;Adjust cosine offsets >"
            OFFSET 1,ChairOff
            WAITC  1,SINEOFF       ;Wait for end of chair cycle >"
            RATE   0,0             ;Stop both cosines  >"
            RATE   1,0             ;                   >"
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOVI   ChairTmp,0      ;Set chair amplitude >"
            JUMP   IDLELOOP        ;Return to idle loop


;-----------------------------------------------------------------------------
; GAP: No stimuli for Gap1Dur+Gap3Dur ms and FlashDur ms light pulse halfway
;-----------------------------------------------------------------------------
GAP:    'G  MOVI   BlockFlg,1      ;Start Gap block    >GAP BLOCK
            DIGOUT [.......0]      ;Reset to initial state >"
            RATE   0,0             ;Stop cosine on drum >"
            RATE   1,0             ;Stop cosine on chair >"
            DAC    0,DrumOff       ;Stop the drum      >"
            DAC    1,ChairOff      ;Stop the chair     >"
            MOV    GapCtr,Gap1Dur  ;Set duration of first half >"
GAP1:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   GapCtr,GAP1     ;Repeat until counter hits zero >"
            BEQ    FlashDur,0,GSKIP ;Skip to GAP2 if FlashDur is zero >"
            DIGOUT [.......1]      ;Turn light on      >"
            MOV    GapCtr,FlashDur ;Set duration of light pulse >"
GAP2:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   GapCtr,GAP2     ;Repeat until counter hits zero >"
GSKIP:      DIGOUT [.......0]      ;Turn light off     >"
            MOV    GapCtr,Gap3Dur  ;Set duration of second half >"
GAP3:       DAC    0,DrumOff       ;Apply drum and chair drift correction >"
            DAC    1,ChairOff      ;                   >"
            DBNZ   GapCtr,GAP3     ;Repeat until counter hits zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; TEST: pre/post test VORD (chair only) block
;-----------------------------------------------------------------------------
TEST:   'P  MOVI   BlockFlg,1      ;Start TEST block   >TESTING
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            MOV    ChairCtr,NchairP ;Set number of cycles to run >"
            DIGOUT [.......0]      ;Ensure light is off >"
            SZ     1,ChairAmp      ;Start chair cosine >"
            OFFSET 1,ChairOff      ;Set cosine offset  >"
            PHASE  1,-90           ;Set cosine relative phase >"
            ANGLE  1,0             ;Set cosine angle   >"
            RATE   1,ChairFrq      ;Set cosine frequency >"
TEST1:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,TEST1         ;Wait for 0 phase   >"
            DBNZ   ChairCtr,TEST1  ;Run cycles until counter hits zero >"
            CLRC   1               ;Stop chair sine at 0 phase >"
TEST2:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,TEST2         ;Wait for end of cycle >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; OPTO SINE Training Block: Turn sinusoidal chair + opto stimulus on/off
;
; APPROACH 1: DELAY from previous cycle (no phase shift needed)
;   - CyclCtr = EveryNth - 1 (set by .s2s), fires from previous cycle's WAITC
;   - Chair phase is standard (-90), angle is 0 (NO phase shift)
;   - DelVal1 = pre-computed adjusted delay (from table, includes period offset)
;   - DelCh1 = channel for first hemisphere (3 or 4, from table)
;   - DelVal2 = halfCycle - 4 (inter-hemisphere overhead, set by .s2s)
;   - Second hemisphere always fires on opposite channel
;
; Instruction path (WAITC fall-through to first DIGPC):
;   DBNZ(1) + BGE(1) + ADDI(1) + BEQ_noStim(1) = 4 pre-DELAY overhead
;   DELAY(DelVal1) = adjusted delay
;   BEQ_ch(1) + DIGPS(1) + DIGPS(1) + DIGPS(1) + DIGPC(1) = 5 post-DELAY
;   Total: 9 + DelVal1 steps from WAITC to first DIGPC
;
; First DIGPC to second DIGPC:
;   DELAY(DelVal2) + DIGPS(1) + DIGPS(1) + DIGPS(1) + DIGPC(1) = 4 + DelVal2
;
;-----------------------------------------------------------------------------
OPTOON: 'T  MOVI   BlockFlg,1      ;Start opto stim trials >TRAINING
            MOVI   DrumTmp,0       ;Set drum amplitude >"
            MOV    ChairTmp,ChairAmp ;Set chair amplitude >"
            DIGOUT [.......0]      ;Ensure opto outputs off >"
            SZ     1,ChairAmp      ;Set cosine amplitude >"
            OFFSET 1,ChairOff      ;Set cosine offset  >"
            PHASE  1,-90           ;Standard cosine phase (no shift) >"
            ANGLE  1,0             ;Standard cosine angle (no shift) >"
            RATE   1,ChairFrq      ;Set cosine frequency >"

OPTO0:      MOV    CyclCtr,EveryNth ;Set cycle counter (EveryNth-1 from .s2s) >"
            TABLD  DelVal1,[DelIdx] ;Get adjusted delay from table >"
            TABLD  DelCh1,[DelIdx+5000] ;Get channel for first hemisphere >"

OPTO1:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,OPTO1         ;Wait for 0 phase   >"
            DBNZ   CyclCtr,OPTO1   ;Repeat until cycle counter hits zero >"
            BGE    DelIdx,NDelays,OPTOOFF ;End block if done >"
            ADDI   StimNum,1       ;Increment StimNum by 1 >"

            BEQ    DelVal1,-1,OPTO2 ;Skip stim (NoStim) if DelVal1 = -1 >"
            DELAY  DelVal1         ;Wait adjusted delay from prev cycle >"
            BEQ    DelCh1,4,H1CH4  ;Select channel for first hemisphere >"

H1CH3:      DIGPS  2,P,PulInt      ;First hemisphere: channel 2 >"
            DIGPS  2,D,PulDur      ;                   >"
            DIGPS  2,C,PulNCycl    ;                   >"
            DIGPC  2,G             ;Start first hemisphere burst >"
            DELAY  DelVal2         ;Inter-hemisphere delay >"
            DIGPS  3,P,PulInt      ;Second hemisphere: channel 3 >"
            DIGPS  3,D,PulDur      ;                   >"
            DIGPS  3,C,PulNCycl    ;                   >"
            DIGPC  3,G             ;Start second hemisphere burst >"
            JUMP   OPTO2           ;                   >"

H1CH4:      DIGPS  3,P,PulInt      ;First hemisphere: channel 3 >"
            DIGPS  3,D,PulDur      ;                   >"
            DIGPS  3,C,PulNCycl    ;                   >"
            DIGPC  3,G             ;Start first hemisphere burst >"
            DELAY  DelVal2         ;Inter-hemisphere delay >"
            DIGPS  2,P,PulInt      ;Second hemisphere: channel 2 >"
            DIGPS  2,D,PulDur      ;                   >"
            DIGPS  2,C,PulNCycl    ;                   >"
            DIGPC  2,G             ;Start second hemisphere burst >"
            JUMP   OPTO2           ;                   >"

OPTO2:      ADDI   DelIdx,1        ;Increment DelIdx by 1 >"
            JUMP   OPTO0

OPTOOFF: 't CLRC   1               ;Stop chair sine at 0 phase >"
OPTO3:      OFFSET 1,ChairOff      ;Adjust cosine offset >"
            WAITC  1,OPTO3         ;Wait for end of cycle >"
            RATE   1,0             ;Stop chair cosine  >"
            MOVI   DrumTmp,0       ;Set drum amplitude to zero >"
            MOVI   ChairTmp,0      ;Set chair amplitude to zero >"
            MOVI   BlockFlg,0      ;Set block as inactive >"
            JUMP   IDLELOOP


;-----------------------------------------------------------------------------
; WAVE: Custom waveform stimuli (if provided)
;-----------------------------------------------------------------------------
;WAVE:   'W  WAVEGO w,TW
;            WAVEST T