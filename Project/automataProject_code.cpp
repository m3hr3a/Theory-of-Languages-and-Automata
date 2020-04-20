#include <bits/stdc++.h>
#include <string>
#include <vector>
#include <map>
#include <stack>

char epsilon = char(-18);
using namespace std;

class CFG{
public:
    string startVar; // start variable
    vector<string> variables; // variables
    vector<string> terminals; // terminals
    vector<string>* grammers; // grammers, for each variable we have a vector of destinations
    int nV, nT; // nV : number of variables , nT : number of terminals

    void checkInput() {
        cout << "Variables :" << endl;
        for (int i = 0; i < nV; i++) {
            cout << variables[i] << endl;
        }
        cout << "terminals : " << endl;
        for (int i = 0; i < nT; i++) {
            cout << terminals[i] << endl;
        }
        cout << "grammers : " << endl;
        for (int i = 0; i < nV; i++) {
            cout << variables[i] << ":" << endl;
            for (int j = 0; j < grammers[i].size(); j++) {
                cout << "-" <<  grammers[i][j] << endl;
            }
        }
    }

    void getInput() {
        string varStr, termStr, temp, temp1, convertTo;
        nV = 0;
        nT = 0;
        cin >> startVar;
        cin.ignore();
        getline(cin, varStr);
        stringstream X1(varStr);
        while (getline(X1, temp, ' ')) {
            nV++;
            variables.push_back(temp);
        }
        getline(cin, termStr);
        stringstream X2(termStr);
        while (getline(X2, temp, ' ')) {
            nT++;
            terminals.push_back(temp);
        }
        grammers = new vector<string> [nV];
        for (int i = 0; i < nV; i++) {
            getline(cin, temp);
            temp = temp.substr(2, temp.length() - 2);
            stringstream X(temp);
            while (getline(X, temp1, '|')) {
                grammers[i].push_back(temp1);
            }
        }
    }

    void toPDA() {
        cout << "qStart" << endl;  // start state
        cout << "" << endl;  // stack is empty at first
        cout << "qAccpet" << endl; // accept state
        for (auto i : terminals) {
            cout << i << " ";
        } cout << endl;
        // alphabet
        // stack alphabet
        for (int i = 0; i < variables.size(); i++) {
            cout << variables[i] << " ";
        }
        for (auto i : terminals) {
            cout << i << " ";
        } cout << endl;
        int c = 0;
        for(int i = 0; i < nV; i++) {
            c += grammers[i].size();
        }
        c += terminals.size();
        c += 2;
        // print states
        cout << "qStart " << "q "  << "qAccept" << endl;
        // print number of ruless
        cout << c << endl;
        cout << "qStart," <<  epsilon << "," << epsilon << ",q," << startVar <<"$" << endl;
        // for each A -> w : e,A -> w
        for(int i = 0; i < nV; i++) {
            for (auto j : grammers[i]) {
                cout << "q," <<  epsilon << "," << variables[i] << ",q," << j << endl;
            }
        }
        // for each a : a,a -> e
        for (auto i : terminals) {
            cout << "q," <<  i << "," << i << ",q," << epsilon << endl;
        }
        cout << "q," <<  epsilon << "," << "$" << ",qAccept," << epsilon << endl;
    }
};

class NFA {
public:

    string startState; // name of NFA's start state
    vector<string> FinalStates; // name of final states
    vector<string> alphabet; // alpahbets
    vector<string> states; // name of states
    map<string, int> state2num; // we store integer number related to each state name here
    // 2 next used for NFA to RegEx
    map<int, string>* delta; // delta[i] stores (j,exp) which there is an edge from i to j by exp
    map<int, string>* deltaInv; // deltaInv[i] stores (j,exp) which there is an edge from j to i by exp
    // next used for NFA to minDFA
    multimap<string, int>* deltaTo; // deltaTo[i] stores all (s,j) : there is an edge from i to j by s
    int nF = 0, nS = 0, nSt = 0; // nF : number of Final States, nS : size of alphabet, nSt : number of states

    void printNFA() {
        //cout << nSt << " " << nS << " " << nF << endl;
        cout << startState << endl;
        for (int i = 0; i < nF; i++) {
            cout << FinalStates[i] << " ";
        }
        cout << endl;
        for (int i = 0; i < nS - 1; i++) {
            cout << alphabet[i] << " ";
        }
        cout << endl;
        for (int i = 1; i < nSt - 1; i++) {
            cout << states[i] << " ";
        }
        cout << endl;
        bool flag = false;
        map<int, string>::iterator itr;
        multimap<string, int>::iterator itr3;
        for (int i = 0; i < nSt - 2; i++) {
            for (auto k : alphabet) {
                auto res = deltaTo[i].equal_range(k);
                int c = 0;
                for (auto itr = res.first; itr!=res.second; itr++) {
                    if (c == 0)
                        cout << states[itr->second+1];
                    else {
                        cout << "," << states[itr->second+1];
                    }
                    c++;
                    flag = true;
                }
                if (!flag) {
                    cout << "-";
                }
                cout << " ";
                flag = false;
            }
            cout << endl;
        }
        for (auto i:state2num) {
            //cout << i.first << " " << i.second << endl;
        }
    }

    void getInput() {
        string temp;
        // get startState
        cin >> startState;
        cin.ignore();
        getline(cin, temp);
        stringstream X1(temp);
        // get final states
        while (getline(X1, temp, ' ')) {
            nF++;
            FinalStates.push_back(temp);
        }
        getline(cin, temp);
        stringstream X2(temp);
        // get alphabet
        while (getline(X2, temp, ' ')) {
            nS++;
            alphabet.push_back(temp);
        }
        alphabet.push_back(""); nS++;
        getline(cin, temp);
        stringstream X3(temp);
        // added start
        state2num.insert(make_pair("qStart", nSt));
        states.push_back("qStart"); nSt++;
        // get states and update numbers related to each state in state2num
        while (getline(X3, temp, ' ')) {
            state2num.insert(make_pair(temp, nSt));
            nSt++;
            states.push_back(temp);
        }
        // added end
        state2num.insert(make_pair("qAccept", nSt));
        states.push_back("qAccept"); nSt++;
        delta = new map <int, string> [nSt];
        deltaInv = new map <int, string> [nSt];
        deltaTo = new multimap<string, int> [nSt-2]; // for DFA usage we have 2 less states!

        // connect new start to NFA start
        delta[0].insert(make_pair(state2num.find(startState)->second, alphabet[nS-1]));
        deltaInv[state2num.find(startState)->second].insert(make_pair(0, alphabet[nS-1]));
        // get transitions
        string temp2;
        for (int i = 1; i < nSt - 1; i++) {
            getline(cin, temp);
            stringstream X4(temp);
            int j = 0;
            while (getline(X4, temp2, ' ')) {
                stringstream X5(temp2);
                while (getline(X5, temp, ',')) {
                    if (temp != "-") {
                        deltaTo[i-1].insert(make_pair(alphabet[j], state2num.find(temp)->second - 1));
                        if (delta[i].find((state2num.find(temp)->second)) == delta[i].end()) {
                            delta[i].insert(make_pair(state2num.find(temp)->second, alphabet[j]));
                        } else {
                            delta[i].find((state2num.find(temp)->second))->second = "(" +
                                                                                    delta[i].find((state2num.find(
                                                                                            temp)->second))->second
                                                                                    + "|" + alphabet[j] + ")";
                        }
                        if (deltaInv[state2num.find(temp)->second].find(i) ==
                            deltaInv[state2num.find(temp)->second].end()) {
                            deltaInv[state2num.find(temp)->second].insert(make_pair(i, alphabet[j]));
                        } else {
                            deltaInv[state2num.find(temp)->second].find(i)->second =
                                    "(" + deltaInv[state2num.find(temp)->second].find(i)->second
                                    + "|" + alphabet[j] + ")";
                        }
                    }
                }
                j++;
            }
        }
        // connect NFAs Finals to added end
        for (int i = 0; i < nF; i++) {
            deltaInv[nSt-1].insert(make_pair(state2num.find(FinalStates[i])->second, alphabet[nS-1]));
            delta[state2num.find(FinalStates[i])->second].insert(make_pair(nSt-1, alphabet[nS-1]));
        }
    }

    void deleteState(int number) {
        string temp;
        multimap <int, string> :: iterator itr1;
        multimap <int, string> :: iterator itr2, itr3;
        string temp1 = "";
        // * for self loop
        for (itr1 = delta[number].begin(); itr1 != delta[number].end(); itr1++) {
            if (itr1->first == number) {
                if (itr1->second.empty()) {
                    temp1 = "(" + string(1, epsilon) + ")*";
                } else {
                    temp1 = "(" + itr1->second + ")*";
                }

            }
        }
        bool ok = false;
        for (itr1 = deltaInv[number].begin(); itr1 != deltaInv[number].end(); itr1++) {
            for (itr2 = delta[number].begin(); itr2 != delta[number].end(); itr2++) {
                if (number != itr1->first & number != itr2->first) {
                    if (states[itr1->first] != "NULL" && states[itr2->first] != "NULL") {
                        temp = itr1->second + temp1 + itr2->second; // concat for new connecting regEx
                        for (itr3 = delta[itr1->first].begin(); itr3 != delta[itr1->first].end(); itr3++) {
                            // in this loop edges will be updated by union for delta
                            if (itr3->first == itr2->first) {
                                if (itr3->second != temp) {
                                    if (temp.empty()) {
                                        itr3->second = "(" + itr3->second + "|" + string(1, epsilon) + ")";
                                    } else if (itr3->second.empty()) {
                                        itr3->second = "(" + string(1, epsilon) + "|" + temp + ")";
                                    } else {
                                        itr3->second = "(" + itr3->second + "|" + temp + ")";
                                    }
                                }
                                ok = true;
                                break;
                            }
                        }
                        if (!ok) {
                            delta[itr1->first].insert(make_pair(itr2->first, temp));
                        }
                        ok = false;
                        for (itr3 = deltaInv[itr2->first].begin(); itr3 != deltaInv[itr2->first].end(); itr3++) {
                            // another step of edge updating for deltaINV
                            if (itr3->first == itr1->first) {
                                if (itr3->second != temp) {
                                    if (temp.empty()) {
                                        itr3->second = "(" + itr3->second + "|" + string(1, epsilon) + ")";
                                    } else if (itr3->second.empty()) {
                                        itr3->second = "(" + string(1, epsilon) + "|" + temp + ")";
                                    } else {
                                        itr3->second = "(" + itr3->second + "|" + temp + ")";
                                    }
                                }
                                ok = true;
                                break;
                            }
                        }
                        if (!ok) {
                            deltaInv[itr2->first].insert(make_pair(itr1->first, temp));
                        }
                        ok = false;
                    }
                }
            }
        }
        states[number] = "NULL";
    }

    void toRegEX() {
        for (int i = 1; i < nSt - 1; i++) {
            deleteState(i);
        }
        cout << delta[0].find(nSt-1)->second;
    }

    void toDFA() {
        map<string, int> dfaStoN;
        int n = nSt-2;
        int nDFA = pow(2,n);
        map<int,int> E[n];
        vector<string> str_DFA_states;
        map<string,int>::iterator itr;
        bool flag = false;

        /** Add each state to E[state] **/
        for (int i = 0; i < n; i++) {
            E[i].insert({i, i});
        }
        /** Calulate E[state]**/
        for (int i = 0; i < n; i++) {
            flag = false;
            while (true) {
                for (auto itr1 = E[i].begin(); itr1 != E[i].end(); itr1++) {
                    flag = false;
                    auto result = deltaTo[itr1->second].equal_range("");
                    for (itr = result.first; itr != result.second; itr++) {
                        if (E[i].find(itr->second) == E[i].end()) {
                            E[i].insert({itr->second, itr->second}); flag = true;
                        }
                    }
                }
                if (!flag) {
                    break;
                }
            }

        }
        /** DFA states **/
        for (int i = 0; i < nDFA; i++) {
            string str = "";
            for (int j = 0; j < n; j++) {
                if ((i & (1 << j)) > 0)
                    if (str == "") {
                        str = to_string(j);
                    }
                    else {
                        str = str + "," + to_string(j);
                    }
            }
            str_DFA_states.push_back(str);
            dfaStoN.insert({str,i});
        }

        /** cal dfa start = E[q0] **/
        auto tempE = E[state2num[startState] - 1];
        vector <int> t2;
        for (auto i = tempE.begin(); i != tempE.end(); i++) {
            t2.push_back(i->second);
        }
        sort(t2.begin(), t2.end());
        string ts  = to_string(t2[0]);
        for (int i = 1; i < t2.size(); i++) {
            ts.append(",").append(to_string(t2[i]));
        }
        int dfaStartNum =  dfaStoN.find(ts)->second;

        vector <int> tvec;
        vector <int> dfaFinal;
        map <int, int> nfaF;
        for (int i = 0; i < FinalStates.size(); i++) {
            nfaF.insert({state2num[FinalStates[i]]-1,state2num[FinalStates[i]]-1});
        }

        string temp;
        int dfaDelta[nDFA][nS-1];

        /** cal dfa's transitions **/
        for (int i = 0; i < nDFA; i++) {
            for (int j = 0; j < nS - 1; j++) {
                tvec.clear();
                flag = false;
                stringstream X(str_DFA_states[i]);
                while (getline(X, temp, ',')) {
                    if (nfaF.find(stoi(temp)) != nfaF.end()) {
                        bool f = true;
                        for (int u = 0; u < dfaFinal.size(); u++) {
                            if (dfaFinal[u] == i) {
                                f = false;
                            }
                        }
                        if (f) {
                            dfaFinal.push_back(i);
                        } f = true;
                    }
                    auto res = deltaTo[stoi(temp)].equal_range(alphabet[j]);
                    for (itr = res.first; itr != res.second; itr++) {
                        for (auto k = E[itr->second].begin(); k != E[itr->second].end(); k++) {
                            for (int o = 0; o < tvec.size(); o++) {
                                if (tvec[o] == k->second) {
                                    flag = true;
                                }
                            }
                            if (!flag) {
                                tvec.push_back(k->second);
                            }
                            flag = false;
                        }
                    }
                }
                sort(tvec.begin(), tvec.end());
                string ts = "";
                for (int l = 0; l < tvec.size(); l++) {
                    if (ts == "")
                        ts = to_string(tvec[0]);
                    else
                        ts.append(",").append(to_string(tvec[l]));
                }
                dfaDelta[i][j] = dfaStoN.find(ts)->second;
            }
        }

        /** delte not reached states **/
        map <int,int> deleted;
        map<int,int> notDel;
        int c2 = 0;
        notDel.insert({dfaStartNum,dfaStartNum});
        while (true) {
            flag = true;
            for (auto itr = notDel.begin(); itr != notDel.end();itr++) {
                for (int k = 0; k < nS-1;k++) {
                    auto s = dfaDelta[itr->second][k];
                    if (notDel.find(s) == notDel.end()) {
                        flag = false;
                        notDel.insert({s,s});
                    }
                }
            }
            if (flag) {
                break;
            }
        }
        for (int i = 0; i < nDFA; i++) {
            if (notDel.find(i) == notDel.end()) {
                if (deleted.find(i) == deleted.end()) {
                    deleted.insert({i,i});
                }
            }
        }

        int eqv[nDFA][nDFA]; // if nDFA gets large, memory error may occor;
        /** minimization alg. **/
        for (int i = 0; i < nDFA; i++) {
            for (int j = i + 1; j < nDFA; j++) {
                eqv[i][j] = 0;
            }
        }
        int c  =0;
        for (int i = 0; i < nDFA; i++) {
            for (int j = i + 1; j < nDFA ;j++) {
                if (deleted.find(i) == deleted.end() && deleted.find(j) == deleted.end()) {
                    bool f1 = false;
                    bool f2 = false;
                    for (int k = 0; k < dfaFinal.size(); k++) {
                        if (dfaFinal[k] == i) {
                            f1 = true;
                        }
                    }
                    for (int k = 0; k < dfaFinal.size(); k++) {
                        if (dfaFinal[k] == j) {
                            f2 = true;
                        }
                    }
                    if (f1 ^ f2) {
                        eqv[i][j] = 1;
                    }
                } else {
                    eqv[i][j] = -1;
                }
            }
        }
        while (true) {
            flag = true;
            for (int i = 0; i < nDFA; i++) {
                for (int j = i + 1; j < nDFA; j++) {
                    if (deleted.find(i) == deleted.end() && deleted.find(j) == deleted.end()){
                        for (int k = 0; k < nS - 1; k++) {
                            if (dfaDelta[i][k] < dfaDelta[j][k]) {
                                if (eqv[dfaDelta[i][k]][dfaDelta[j][k]] > 0) {
                                    if (eqv[i][j] == 0) {
                                        eqv[i][j] = 1;
                                        flag = false;
                                    }
                                }
                            } else if (dfaDelta[i][k] > dfaDelta[j][k]){
                                if (eqv[dfaDelta[j][k]][dfaDelta[i][k]] > 0) {
                                    if (eqv[i][j] == 0) {
                                        eqv[i][j] = 1;
                                        flag = false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (flag) {
                break;
            }
        }

        /** cal min. dfa parameters **/
        int group[nDFA];
        for (int i = 0; i < nDFA; i++) {
            group[i] = -1;
        }

        c = 0;

        for (int i = 0; i < nDFA; i++) {
            if (group[i] == -1 && deleted.find(i) == deleted.end()) {
                group[i] = c; c++;
            }
            for (int j = i + 1; j < nDFA ;j++) {
                if (eqv[i][j] == 0) {
                    group[j] = c - 1;
                }
            }
        }

        int minDFA_start = group[dfaStartNum];
        map <int, int> minDFA_final;
        for (int i : dfaFinal) {
            if (group[i] != -1) {
                minDFA_final.insert({group[i], group[i]});
            }
        }

        int minDFA_delta[c][nS - 1];
        for (int i = 0; i < c; i++) {
            for (int j = 0;  j < nS - 1; j++) {
                for (int k = 0; k < nDFA; k++) {
                    if (group[k] == i) {
                        minDFA_delta[i][j] = group[dfaDelta[k][j]];
                        break;
                    }
                }
            }
        }

        /** print output :: minimized DFA  **/
        cout << "q" + to_string(minDFA_start) << endl;

        for (auto & i : minDFA_final) {
            cout << "q" + to_string(i.second) + " ";
        } cout << endl;

        for (int i = 0; i < nS-1;i++) {
            cout << alphabet[i] << " ";
        } cout << endl;

        for (int i  = 0; i < c; i++) {
            cout << "q" + to_string(i) + " ";
        } cout << endl;

        for (int i = 0; i < c; i++) {
            for (int j = 0;  j < nS - 1; j++) {
                cout << "q" + to_string(minDFA_delta[i][j]) << " ";
            }
            cout << endl;
        }
    }

};

NFA* toNFA(string regEx) {
    stack <char> operators;
    stack <int> s1, s2;
    map <int,int> NFAstates;
    map <char,char> alphabets;
    char delta[100][100];
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < 100; j++) {
            delta[i][j] = '*';
        }
    }
    int nS = 0;
    for (int i = 0; i < regEx.length(); i++) {
        // union sign
        if (regEx[i] == '|') {
            nS++;
            operators.push(regEx[i]);
        }
        // (
        if (regEx[i] == '(') {
            nS++;
            operators.push('(');
            s1.push(-1);
            s2.push(-1);
        }
        // )
        if (regEx[i] == ')') {
            int start, end; bool u = false;
            // if union build new union block
            if (operators.top() == '|') {
                u = true;
                start = ++nS;  end = ++nS;
                NFAstates.insert({start,start});
                NFAstates.insert({s1.top(),s1.top()});
                NFAstates.insert({s2.top(),s2.top()});
                NFAstates.insert({end,end});
                delta[start][s1.top()] = epsilon;
                delta[s2.top()][end] = epsilon;
                s1.pop(); s2.pop();
            } // add more to new union block
            while (operators.top() == '|') {
                operators.pop();
                NFAstates.insert({start,start});
                NFAstates.insert({s1.top(),s1.top()});
                NFAstates.insert({s2.top(),s2.top()});
                NFAstates.insert({end,end});
                delta[start][s1.top()] = epsilon;
                delta[s2.top()][end] = epsilon;
                s1.pop(); s2.pop();
            }
            // add new built block
            if (u) {
                s1.push(start); s2.push(end);
            }

            if (operators.top() == '(') {
                // remove old states
                int endstate = s2.top();
                int t = s1.top(); s2.pop(); s1.pop();
                while (s1.top() != -1) {
                    NFAstates.insert({s2.top(),s2.top()});
                    NFAstates.insert({t,t});
                    delta[s2.top()][t] = epsilon;
                    t = s1.top(); s2.pop(); s1.pop();
                }
                if (s1.top() == -1) {
                    s1.pop(); s2.pop();}

                s1.push(t); s2.push(endstate);
                if (i < regEx[i+1] - 1 && regEx[i+1] == '*') {
                    int newStart = ++nS;
                    int newEnd = ++nS;
                    NFAstates.insert({newStart,newStart});
                    NFAstates.insert({s1.top(),s1.top()});
                    NFAstates.insert({s2.top(),s2.top()});
                    NFAstates.insert({newEnd,newEnd});
                    delta[s2.top()][s1.top()] = epsilon;
                    delta[newStart][s1.top()] = epsilon;
                    delta[s2.top()][newEnd] = epsilon;
                    delta[newStart][newEnd] = epsilon;
                    s1.pop(); s2.pop();
                    s1.push(newStart);
                    s2.push(newEnd);
                    operators.pop();
                } else {
                    operators.pop();
                }
            }
            u = false;
        }
        // build simple alphabet concat
        if (regEx[i] != ')' && regEx[i] != '*' && regEx[i] != '(' && regEx[i] != '|') {
            nS++;
            //cout << nS << " " << nS+1 << " " << regEx[i] << endl;
            alphabets.insert({regEx[i], regEx[i]});
            NFAstates.insert({nS,nS});
            NFAstates.insert({nS+1,nS+1});
            s1.push(nS);
            delta[nS][nS+1] = regEx[i]; nS++;
            while (i < regEx.length() - 1 && regEx[i+1] != ')' && regEx[i+1] != '*' && regEx[i+1] != '(' && regEx[i+1] != '|') {
                //cout << nS << " " << nS+1 << " " << regEx[i+1] << endl;
                NFAstates.insert({nS,nS});
                alphabets.insert({regEx[i+1], regEx[i+1]});
                NFAstates.insert({nS+1,nS+1});
                delta[nS][nS+1] = regEx[i+1]; nS++; i++;
            }
            s2.push(nS);
        }
    }
    // end
    int endstate = s2.top();
    int t = s1.top(); s2.pop(); s1.pop();
    while (!s1.empty()) {
        delta[s2.top()][t] = epsilon;
        t = s1.top(); s2.pop(); s1.pop();
    }
    // fix output as a NFA
    auto myNFA = new NFA();
    int NFA_start = t;
    int NFA_end = endstate;
    vector <char> DFA_alph;
    myNFA->startState = "q" + to_string(NFA_start);
    myNFA->FinalStates.push_back("q" + to_string(NFA_end));
    for (auto itr : alphabets) {
        myNFA->nS++;
        DFA_alph.push_back(itr.second);
        myNFA->alphabet.push_back(string(1,itr.second));
    }
    myNFA->alphabet.push_back(""); myNFA->nS++;
    DFA_alph.push_back(epsilon);
    myNFA->states.push_back("qStart"); myNFA->nSt++;
    for (auto itr : NFAstates) {
        myNFA->nSt++;
        myNFA->state2num.insert({"q" + to_string(itr.second), myNFA->nSt - 1});
        myNFA->states.push_back("q" + to_string(itr.second));
    } myNFA->nSt++;
    myNFA->nF = 1;
    myNFA->deltaTo = new multimap<string,int>[myNFA->nSt-2];
    bool flag = true;
    int c1 = -1, c2 = -1;
    for (auto i : NFAstates) {
        c1++;
        for (auto k : DFA_alph) {
            c2 = -1;
            for (auto j : NFAstates) {
                c2++;
                if (delta[i.second][j.second] == k) {
                    if (k == epsilon) {
                        myNFA->deltaTo[c1].insert({"",c2});
                    }
                    else
                        myNFA->deltaTo[c1].insert({string(1,k),c2});
                }
            }
            if (flag) {
            }
            flag = true;
        }
    }
    return myNFA;
}

int main() {
    string inputMode, outputMode;
    cin >> inputMode;
    string RegEx;

    if (inputMode == "RegEx") {
        cin >> RegEx;
        auto myNFA = toNFA(RegEx);
        cin >> outputMode;
        if (outputMode == "NFA")  {
            myNFA->printNFA();
        }
        if (outputMode == "DFA")  {
            myNFA->toDFA();
        }
    }

    if (inputMode == "NFA") {
        auto myNfa = new NFA();
        myNfa->getInput();
        cin >> outputMode;
        if (outputMode == "RegEx")  {
            myNfa->toRegEX();
        }
        if (outputMode == "DFA")  {
            myNfa->toDFA();
        }
    }

    if (inputMode == "CFG") {
        auto myCFG = new CFG;
        myCFG->getInput();
        cin >> outputMode;
        if (outputMode == "PDA") {
            myCFG->toPDA();
        }
    }

}
