#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "iostream"

void drawHistBox(double x1=0., double x2=1., double y1=0., double y2=1., int lineWidth=3333, int lineStyle=1, int lineColor=1)
//void drawHistBox(double x1, double x2, double y1, double y2, int lineWidth)
{
  TLine *l1 = new TLine(x1,y1,x2,y1);
  l1->SetLineWidth(lineWidth/1000);
  l1->SetLineStyle(lineStyle);
  l1->SetLineColor(lineColor);
  l1->Draw("same");
  TLine *l2 = new TLine(x1,y2,x2,y2);
  l2->SetLineWidth((lineWidth%1000)/100);
  l2->SetLineStyle(lineStyle);
  l2->SetLineColor(lineColor);
  l2->Draw("same");
  TLine *l3 = new TLine(x1,y1,x1,y2);
  l3->SetLineWidth((lineWidth%100)/10);
  l3->SetLineStyle(lineStyle);
  l3->SetLineColor(lineColor);
  l3->Draw("same");
  TLine *l4 = new TLine(x2,y1,x2,y2);
  l4->SetLineWidth(lineWidth%10);
  l4->SetLineStyle(lineStyle);
  l4->SetLineColor(lineColor);
  l4->Draw("same");
}

void drawText(double x1=0., double y1=0., const char* text="test", int textFont=42, double textSize=0.05, double textAngle=0, int textColor=1)
//void drawText(double x1, double y1, const char* text, int textFont, double textSize, double textAngle)
{
  TLatex *tex = new TLatex(x1, y1, text);
  tex->SetTextFont(textFont);
  tex->SetTextSize(textSize);
  tex->SetTextAngle(textAngle);
  tex->SetTextColor(textColor);
  tex->Draw("same");
}

void drawLine(double x1=0., double y1=0., double x2=0., double y2=1., int lineWidth=1, int lineStyle=1, int lineColor=1)
//void drawLine(double x1, double y1, double x2, double y2, int lineWidth, int lineStyle, int lineColor)
{
  TLine *la = new TLine(x1,y1,x2,y2);
  la->SetLineWidth(lineWidth);
  la->SetLineStyle(lineStyle);
  la->SetLineColor(lineColor);
  la->Draw("same");
}

void drawSysError(TGraphErrors *gr, double xoffset=0.05, double yoffset=0.03, int lineColor=1, bool logx=0, bool logy=0)
{
  if(!gr) {
    std::cout << "No TGraphErrors, return!" << std::endl;
    return;
  }
  
  for(int i=0;i<gr->GetN();i++) {
    double x0 = gr->GetX()[i];
    double y0 = gr->GetY()[i];
    double ye = gr->GetEY()[i];
    double x1 = x0-xoffset;
    double x2 = x0+xoffset;
    double y1 = y0-ye;
    double y2 = y0+ye;
    double y11 = y1+yoffset;
    double y22 = y2-yoffset;
    if(logx) {
      x1 = x0*(1-xoffset);
      x2 = x0*(1+xoffset);
    }
    if(logy) {
      y11 = y1*(1+yoffset);
      y22 = y2*(1-yoffset);
    }
    
    TLine *la = new TLine(x1,y1,x2,y1);
    la->SetLineWidth(2);
    la->SetLineColor(lineColor);
    la->Draw("same");
    TLine *la1 = new TLine(x1,y1,x1,y11);
    la1->SetLineWidth(1);
    la1->SetLineColor(lineColor);
    la1->Draw("same");
    TLine *la2 = new TLine(x2,y1,x2,y11);
    la2->SetLineWidth(1);
    la2->SetLineColor(lineColor);
    la2->Draw("same");
    
    TLine *lb = new TLine(x1,y2,x2,y2);
    lb->SetLineWidth(2);
    lb->SetLineColor(lineColor);
    lb->Draw("same");
    TLine *lb1 = new TLine(x1,y2,x1,y22);
    lb1->SetLineWidth(1);
    lb1->SetLineColor(lineColor);
    lb1->Draw("same");
    TLine *lb2 = new TLine(x2,y2,x2,y22);
    lb2->SetLineWidth(1);
    lb2->SetLineColor(lineColor);
    lb2->Draw("same");
  }
}
