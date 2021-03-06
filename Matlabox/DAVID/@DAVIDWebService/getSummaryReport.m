function xReturn = getSummaryReport(obj)
%getSummaryReport(obj)
%
%     Input:
%   
%     Output:
%       return{:} = ()

% Build up the argument lists.
values = { ...
   };
names = { ...
   };
types = { ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://service.session.sample', ...
    'getSummaryReport', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getSummaryReport', ...
    soapMessage);
xReturn = parseSoapResponse(response);
