function xReturn = getDefaultCategoryNames(obj)
%getDefaultCategoryNames(obj)
%
%     Input:
%   
%     Output:
%       return = (string)

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
    'getDefaultCategoryNames', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'urn:getDefaultCategoryNames', ...
    soapMessage);
xReturn = parseSoapResponse(response);
